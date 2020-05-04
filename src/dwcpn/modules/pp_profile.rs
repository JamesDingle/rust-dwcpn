use crate::dwcpn::modules::config::{DEPTH_PROFILE_COUNT, NUM_WAVELENGTHS, WAVELENGTHS, AW, DEPTH_PROFILE_STEP};
use crate::dwcpn::modules::absorption::calc_ac;


pub struct PpProfile {
    pub pp_profile: [f64; DEPTH_PROFILE_COUNT],
    pub par_profile: [f64; DEPTH_PROFILE_COUNT],
    pub euphotic_depth: f64,
    pub euphotic_depth_index: usize,
    pub success: bool;
}

pub fn calculate_bw() -> [f64; NUM_WAVELENGTHS] {

    // scattering coefficient of pure seawater at 500nm
    const BW500: f64 = 0.00288;

    let mut bw: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];

    for i in 0..NUM_WAVELENGTHS {
        bw[i] = BW500 * (WAVELENGTHS[i] / 500.0).powf(-4.3);
    }

    bw
}

pub fn calculate_bbr() -> [f64; NUM_WAVELENGTHS] {
    const BR488: f64 = 0.00027;

    let mut bbr: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];

    for i in 0..NUM_WAVELENGTHS {
        bbr[i] = 0.5 * BR488 * (WAVELENGTHS[i] / 488.0).powf(-5.3);
    }

    bbr
}

pub fn calculate_ay() -> [f64; NUM_WAVELENGTHS] {
    let mut ay: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];

    for i in 0..NUM_WAVELENGTHS {
        ay[i] = (-0.014 * (WAVELENGTHS[i] - 440.0)).exp();
    }

    ay
}

pub fn compute_pp_depth_profile(
    chl_profile: [f64; DEPTH_PROFILE_COUNT],
    depth_profile: [f64; DEPTH_PROFILE_COUNT],
    zenith_r: f64,
    direct_irradiance: [f64; NUM_WAVELENGTHS],
    diffuse_irradiance: [f64; NUM_WAVELENGTHS],
    bw: [f64; NUM_WAVELENGTHS],
    bbr: [f64; NUM_WAVELENGTHS],
    ay: [f64; NUM_WAVELENGTHS],
    province_alpha: f64,
    province_pmb: f64,
    yellow_substance: f64
) -> PpProfile {
    let mut pp_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];
    let mut par_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];
    let mut euphotic_depth: f64 = 0.0;
    let mut euphotic_depth_index: usize = 0;
    let mut success = false;



    let mut i_zero: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];
    let mut mu_d: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];
    let mut i_z: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];
    let mut k: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];

    let zenith_w: f64 = (zenith_r.sin() / 1.333).asin();

    for l in 0..NUM_WAVELENGTHS {
        i_zero[l] = direct_irradiance[l] + diffuse_irradiance[l];
        mu_d[l] = (direct_irradiance[l] * zenith_w.cos() + diffuse_irradiance[l] * 0.831000) / i_zero[l];
        i_z[l] = i_zero[l];
    }

    for z in 0..DEPTH_PROFILE_COUNT {
        let chl = chl_profile[z];

        // TODO: potentially rename to ac_mean
        let (ac, al_mean) = calc_ac(chl);

        let ac440 = ac[8]; // TODO: fix this to search for/interpolate to 440nm

        let power = -chl.log10();
        let ay440 = yellow_substance * ac440;

        let bc660 = 0.407 * chl.powf(0.795);
        let mut bbtilda = (0.78 + 0.42 * power) * 0.01;

        if bbtilda < 0.0005 {
            bbtilda = 0.0005;
        } else if bbtilda > 0.01 {
            bbtilda = 0.01;
        }

        let mut alpha_b = 0.0;

        for l in 0..NUM_WAVELENGTHS {
            wl = WAVELENGTHS[l];
            let a = AW[l] + ac[l] + ay440 * ay[l] + 2.0 * bbr[l];
            let mut bc = bc660 * (660.0 / wl).powf(power);

            if bc < 0.0 { bc = 0.0 };

            let bb = bc * bbtilda + bw[l] * 0.50;

            k[l] = (a + bb) / mu_d[l];

            par_profile[z] = par_profile[z] + i_z[l] * 5.0;

            let x = province_alpha * ac[l] * 6022.0 / (2.77 * 36.0 * al_mean);

            alpha_b = alpha_b + x * 5.0 * i_z[l] / mu_d[l];
            i_z[l] = i_z[l] * (-k[l] * DEPTH_PROFILE_STEP).exp();
        }

        pp_profile[z] = (alpha_b / (1.0 + (alpha_b / province_pmb).powf(2.0)).sqrt()) * chl;

        if z > 0 {
            if par_profile[z] < (0.01 * par_profile[0]) {
                euphotic_depth_index = z - 1;
                euphotic_depth = depth_profile[z - 1] + DEPTH_PROFILE_STEP * (100 * par_profile[z - 1] / par_profile[0]).ln() / (par_profile[z - 1] / par_profile[z]).ln();
                success = true
            }
        }

    }


    return PpProfile{
        pp_profile,
        par_profile,
        euphotic_depth,
        euphotic_depth_index,
        success
    }
}