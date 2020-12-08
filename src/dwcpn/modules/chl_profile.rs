use crate::dwcpn::modules::config::{DEPTH_PROFILE_COUNT, DEPTH_PROFILE_START, DEPTH_PROFILE_STEP};
use std::f64::consts::PI;

const TAU: f64 = PI * 2.0;

pub fn gen_chl_profile(
    surface_chl: f64,
    sigma: f64,
    rho: f64,
    z_m: f64,
) -> ([f64; DEPTH_PROFILE_COUNT], [f64; DEPTH_PROFILE_COUNT]) {
    let b_0: f64 = surface_chl
        / (1.0 + (rho / (1.0 - rho)) * (-z_m.powf(2.0) / (2.0 * sigma.powf(2.0))).exp());

    // research if we will always precompute h or keep the calculation in here
    let h: f64 = sigma * (rho / (1.0 - rho)) * b_0 * TAU.sqrt();
    let gauss_height = h / (sigma * TAU.sqrt());

    let depth_array: [f64; DEPTH_PROFILE_COUNT] = gen_depth_array();

    let mut chl_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];

    for i in 0..DEPTH_PROFILE_COUNT {
        let c = -0.5 * ((depth_array[i] - z_m) / sigma).powf(2.0);
        if c.abs() <= 675.0 {
            chl_profile[i] = gauss_height * c.exp() + b_0;
        } else {
            chl_profile[i] = 0.0;
        }
    }

    return (depth_array, chl_profile);
}

fn gen_depth_array() -> [f64; DEPTH_PROFILE_COUNT] {
    let mut depth_array: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];

    for i in 0..DEPTH_PROFILE_COUNT {
        let i_f64 = i as f64;
        depth_array[i] = DEPTH_PROFILE_START + i_f64 * DEPTH_PROFILE_STEP;
    }

    return depth_array;
}
