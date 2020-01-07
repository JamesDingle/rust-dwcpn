use crate::dwcpn::modules::config::{DEPTH_PROFILE_START, DEPTH_PROFILE_COUNT, DEPTH_PROFILE_STEP};
use std::f64::consts::PI;

const TAU: f64 = PI * 2.0;

pub struct ChlProfile {
    pub depth_array: [f64; DEPTH_PROFILE_COUNT],
    pub chl_profile: [f64; DEPTH_PROFILE_COUNT]
}

impl ChlProfile {
    pub fn print_profile(&self) {
        for i in 0..DEPTH_PROFILE_COUNT {
            println!("{:?},{:?}", i, self.chl_profile[i]);
        }
    }
}

pub fn gen_chl_profile(surface_chl: f64, sigma: f64, rho: f64, z_m: f64, h: f64) -> ([f64; DEPTH_PROFILE_COUNT], [f64; DEPTH_PROFILE_COUNT]) {
    let gauss_height = h / ( sigma * TAU.sqrt() );
    let b_0: f64 = surface_chl / (1.0 + (rho / (1.0 - rho)) * ( -z_m.powf(2.0) / (2.0 * sigma.powf(2.0)) ).exp() );
    let h: f64 = sigma * (rho / (1.0 - rho)) * b_0 * TAU.sqrt();
    let depth_array: [f64; DEPTH_PROFILE_COUNT] = gen_depth_array();


    let mut chl_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];

    for i in 0..DEPTH_PROFILE_COUNT {
        let c = -0.5 * ( (depth_array[i] - z_m) / sigma ).powf(2.0);
        if c.abs() <= 675.0 {
            chl_profile[i] = gauss_height * c.exp() + b_0;
        } else {
            chl_profile[i] = 0.0;
        }
    }

    return (depth_array, chl_profile)
}

fn gen_depth_array() -> [f64; DEPTH_PROFILE_COUNT] {
    let mut depth_array: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];

    let count = DEPTH_PROFILE_COUNT as f64;
    let step = DEPTH_PROFILE_STEP as f64;

    for i in 0..DEPTH_PROFILE_COUNT {
        let i_f64 = i as f64;
        depth_array[i] = DEPTH_PROFILE_START + i_f64 * DEPTH_PROFILE_STEP;
    }

    return depth_array;
}
