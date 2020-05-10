use crate::dwcpn::modules::chl_profile::gen_chl_profile;
use crate::dwcpn::modules::time::{compute_sunrise, generate_time_array};
use crate::dwcpn::modules::config::{TIMESTEPS, NUM_WAVELENGTHS, WAVELENGTHS, DEPTH_PROFILE_COUNT, DEPTH_PROFILE_STEP};
use crate::dwcpn::modules::zenith::generate_zenith_array;
use crate::dwcpn::modules::irradiance::{compute_irradiance_components, lookup_thekaekara_correction};
use std::f64::consts::PI;
use crate::dwcpn::modules::pp_profile::compute_pp_depth_profile;

pub struct InputParams {
    pub lat: f64,
    pub lon: f64,
    pub z_bottom: f64,
    pub iday: u16,
    pub alpha_b: f64,
    pub pmb: f64,
    pub z_m: f64,
    pub chl: f64,
    pub rho: f64,
    pub h: f64,
    pub sigma: f64,
    pub cloud: f64,
    pub yel_sub: f64,
    pub par: f64,
    pub bw: [f64; NUM_WAVELENGTHS],
    pub bbr: [f64; NUM_WAVELENGTHS],
    pub ay: [f64; NUM_WAVELENGTHS],
}

pub fn calc_pp(input: InputParams) -> (f64, f64) {
    // generate chl depth profile
    let (depth_array, chl_profile) = gen_chl_profile(input.chl, input.sigma, input.rho, input.z_m, input.h);

    // compute sunrise and generate time array
    let (sunrise, delta, phi) = compute_sunrise(input.iday, input.lat);

    let (time_array, delta_t) = generate_time_array(sunrise);
    let (zenith_array, zenith_d_array) = generate_zenith_array(time_array, delta, phi);



    // loop over time array (from sunrise to noon)
    let mut start_time_idx: f64 = -1.0;
    let mut day_length: f64 = 0.0;
    let mut iom: f64 = 0.0;
    let mut delta_prestart: f64 = 0.0;

    let solar_correction = lookup_thekaekara_correction(input.iday);

    let mut surface_irradiance: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let mut par_surface_irradiance: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let mut adjustment: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let mut i_zero: [f64; TIMESTEPS] = [0.0; TIMESTEPS];

    // arrays to store results
    let mut pp: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let mut euphotic_depth: [f64; TIMESTEPS] = [0.0; TIMESTEPS];


    for t in 0..TIMESTEPS {

        // if the zenith angle is yet to go below 80° then skip to the next time step
        if zenith_d_array[t] > 80. {
            continue
        }

        // update start_time if necessary and assign delta_prestart variable
        // to be used for daily integration purposes later.
        // t_start is t_idx value when Zenith angle becomes >80 degrees
        // start_time is calculation start time, start_time_index is the index
        // delta_prestart is time elapsed between dawn and start_time
        let start_time: f64;
        if start_time_idx < 0.0 {
            start_time_idx = t as f64;

            // this line has start_time_idx corrected with +1 as the original code
            // came from FORTRAN which uses indices that start at 1
            start_time = sunrise + delta_t * (start_time_idx + 1.0);

            day_length = 2.0 * (12.0 - (start_time - delta_t));
            iom = input.par * PI / (2.0 * day_length);
            delta_prestart = start_time - sunrise;

        }

        // compute direct and diffuse irradiance components at sea level
        let (mut direct, mut diffuse) = compute_irradiance_components(zenith_array[t], zenith_d_array[t]);

        let mut direct_integrated: f64 = 0.0;
        let mut diffuse_integrated: f64 = 0.0;

        for l in 0..NUM_WAVELENGTHS {
            // apply fractional correction to diffuse and direct components of irradiance
            // the max correction value is 1353.0, so this converts it to as though we were applying a percentage correction
            direct[l] = direct[l] * solar_correction / 1353.0;
            diffuse[l] = diffuse[l] * solar_correction / 1353.0;

            // add this value to the integrated direct/diffuse components
            direct_integrated = direct_integrated + (direct[l] * zenith_array[t].cos());
            diffuse_integrated = diffuse_integrated + diffuse[l];
        }

        surface_irradiance[t] = surface_irradiance[t] + direct_integrated + diffuse_integrated;

        // cloud effect calculations
        let albedo = 0.28 / (1.0 + 6.43 * zenith_array[t].cos());
        let cc = input.cloud / 100.0;
        let idir1 = direct_integrated * (1.0 - cc);
        let flux = ((1.0 - 0.5 * cc) * (0.82 - albedo * (1.0 - cc)) * zenith_array[t].cos()) / ((0.82 - albedo) * zenith_array[t].cos());
        let idif1 = surface_irradiance[t] * flux - idir1;
        let dir_div = idir1 / direct_integrated;
        let dif_div = idif1 / diffuse_integrated;

        for l in 0..NUM_WAVELENGTHS {
            direct[l] = direct[l] + dir_div;
            diffuse[l] = diffuse[l] + dif_div;
        }

        // calculate reflection and convert watts/micron into einsteins/hr/nm
        let zenith_w = (zenith_array[t].sin() / 1.333).asin();
        let mut reflection = 0.5 * (zenith_array[t] - zenith_w).sin().powi(2) / (zenith_array[t] + zenith_w).sin().powi(2);
        reflection = reflection + 0.5 * (zenith_array[t] - zenith_w).tan().powi(2) / (zenith_array[t] + zenith_w).tan().powi(2);

        // recompute surface irradiance across spectrum
        surface_irradiance[t] = 0.0;

        for l in 0..NUM_WAVELENGTHS {
            let wl_coefficient = WAVELENGTHS[l] * 36.0 / (19.87 * 6.022 * 10e6);
            direct[l] = direct[l] * wl_coefficient * zenith_array[t].cos();
            diffuse[l] = diffuse[l] * wl_coefficient;

            surface_irradiance[t] = surface_irradiance[t] + direct[l] + diffuse[l];

            direct[l] = direct[l] * (1.0 - reflection);
            diffuse[l] = diffuse[l] * 0.945;
        }

        // compute surface irradiance from total daily surface irradiance (e.g. satellite par)
        par_surface_irradiance[t] = iom * (PI * (time_array[t] - sunrise) / day_length).sin();
        surface_irradiance[t] = surface_irradiance[t] * 5.0;

        // Adjustment to the difuse and direct component: from use of measured total daily surface irradiance (
        // e.g. satellite PAR) to compute the surface irradiance at all time. SSP
        adjustment[t] = par_surface_irradiance[t] / surface_irradiance[t];

        //compute the adjusted irradiance surface value
        i_zero[t] = 0.0;
        for l in 0..NUM_WAVELENGTHS {
            direct[l] = direct[l] * adjustment[t];
            diffuse[l] = diffuse[l] * adjustment[t];

            i_zero[t] = i_zero[t] + (direct[l] + diffuse[l]) * 5.0;
        }

        let mut pp_profile = compute_pp_depth_profile(
            chl_profile,
            depth_array,
            zenith_array[t],
            direct,
            diffuse,
            input.bw,
            input.bbr,
            input.ay,
            input.alpha_b,
            input.pmb,
            input.yel_sub
        );

        if pp_profile.success == true {
            // clamp euphotic depth to bathymetry if it was calculated higher (in extreme clear water)
            if pp_profile.euphotic_depth.abs() > input.z_bottom.abs() {
                pp_profile.euphotic_depth_index = DEPTH_PROFILE_COUNT - 1;
                pp_profile.euphotic_depth = input.z_bottom;
            }

            euphotic_depth[t] = pp_profile.euphotic_depth;

            // if pp_profile.euphotic_depth_index == 0 {
            //     pp_profile.euphotic_depth_index == 1;
            // }

            for z in 0..pp_profile.euphotic_depth_index {
                pp[t] = pp[t] + DEPTH_PROFILE_STEP * (pp_profile.pp_profile[z] + pp_profile.pp_profile[z+1]) / 2.0;
            }

            pp[t] = pp[t] + pp_profile.pp_profile[pp_profile.euphotic_depth_index] * (euphotic_depth[t] - (pp_profile.euphotic_depth_index as f64 - 1.0) * DEPTH_PROFILE_STEP);


        }

    }

    let mut pp_day = pp[0] * delta_prestart / 2.0;
    // let mut z_phot_day = euphotic_depth[0] * i_zero[0] * delta_prestart / 2.0;
    // let mut i_zero_day = i_zero[0] * delta_prestart / 2.0;

    let mut max_euphotic_depth = 0.0;

    for t in 0..TIMESTEPS-1 {
        pp_day = pp_day + ((pp[t] + pp[t+1]) * delta_t / 2.0);
        if max_euphotic_depth < euphotic_depth[t] {
            max_euphotic_depth = euphotic_depth[t];
        }
    }

    // mutliply by two because we have only integrated over half of the day
    pp_day = pp_day * 2.0;
    return (pp_day, max_euphotic_depth)
}