use crate::dwcpn::modules::chl_profile::{ChlProfile, gen_chl_profile};
use crate::dwcpn::modules::time::{compute_sunrise, generate_time_array};
use crate::dwcpn::modules::config::TIMESTEPS;
use crate::dwcpn::modules::zenith::generate_zenith_array;

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

    let time_array: [f64; TIMESTEPS] = generate_time_array(sunrise);
    let zenith_array: [f64; TIMESTEPS] = generate_zenith_array(time_array, delta, phi);
    0.0
}