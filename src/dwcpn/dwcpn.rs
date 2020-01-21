use crate::dwcpn::modules::chl_profile::gen_chl_profile;
use crate::dwcpn::modules::time::{compute_sunrise, generate_time_array};
use crate::dwcpn::modules::config::TIMESTEPS;
use crate::dwcpn::modules::zenith::generate_zenith_array;
use std::f64::consts::PI;

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
}

pub fn calc_pp(input: InputParams) -> f64 {
    // generate chl depth profile
    let (depth_array, chl_profile) = gen_chl_profile(input.chl, input.sigma, input.rho, input.z_m, input.h);

    // compute sunrise and generate time array
    let (sunrise, delta, phi) = compute_sunrise(input.iday, input.lat);

    let (time_array, delta_t) = generate_time_array(sunrise);
    let (zenith_array, zenith_d_array) = generate_zenith_array(time_array, delta, phi);

    // loop over time array (from sunrise to noon)
    let mut start_time_idx: f64 = -1.0;

    for t in 0..TIMESTEPS {
        let time = time_array[t];

        // if the zenith angle is yet to go below 80° then skip to the next time step
        if zenith_d_array[t] > 80. {
            continue
        }

        // update start_time if necessary and assign delta_prestart variable
        // to be used for daily integration purposes later.
        // t_start is t_idx value when Zenith angle becomes >80 degrees
        // start_time is calculation start time, start_time_index is the index
        // delta_prestart is time elapsed between dawn and start_time
        if start_time_idx < 0.0 {
            start_time_idx = t as f64;

            // this line has start_time_idx corrected with +1 as the original code
            // came from FORTRAN which uses indices that start at 1
            let start_time: f64 = sunrise + delta_t * (start_time_idx + 1.0);

            let day_length: f64 = 2.0 * (12.0 - (start_time - delta_t));

            let iom: f64 = input.par * PI / (2.0 * day_length);

            let delta_prestart: f64 = start_time - sunrise;

        }

        // compute direct and diffuse irradiance components at sea level


    }

    0.0
}