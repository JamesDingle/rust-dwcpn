//use crate::modules::time::gen_time_array;

use crate::dwcpn::dwcpn::{InputParams, calc_pp};
use crate::dwcpn::modules::pp_profile::{calculate_bw, calculate_bbr, calculate_ay};

mod dwcpn;

fn main() {

    // pre calculate bw/bbr/ay arrays for use in the pp_profile calculation later
    // this only needs to be done once for all pixels so I have it here, before we start
    // looping over pixels
    let bw = calculate_bw();
    let bbr = calculate_bbr();
    let ay = calculate_ay();


    let input = InputParams{
        lat: 31.5,
        lon: -64.0,
        z_bottom: -4430.55615573459,
        iday: 5,
        alpha_b: 0.1893,
        pmb: 4.662,
        z_m: 74.7,
        chl: 1.0,
        rho: 0.1,
        sigma: 1.5,
        cloud: 0.0,
        yel_sub: 0.3,
        par: 20.8093246830834,
        bw,
        bbr,
        ay
    };


    println!("pp = {:?}", calc_pp(input));

}
