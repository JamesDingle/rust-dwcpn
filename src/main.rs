//use crate::modules::time::gen_time_array;

use crate::dwcpn::dwcpn::{InputParams, calc_pp};
use crate::dwcpn::modules::pp_profile::{calculate_bw, calculate_bbr, calculate_ay};

mod dwcpn;

fn main() {

    // pre calculate bw/bbr/ay arrays for use in the pp_profile calculation later
    // this only needs to be done once for all pixels so I have it here, before we start
    // looping over pixels
    bw = calculate_bw();
    bbr = calculate_bbr();
    ay = calculate_ay();


    let input = InputParams{
        lat: 45.0,
        lon: 45.0,
        z_bottom: 100.0,
        iday: 100,
        alpha_b: 0.1,
        pmb: 1.0,
        z_m: 20.0,
        chl: 0.1,
        rho: 0.1,
        h: 0.1,
        sigma: 0.1,
        cloud: 0.0,
        yel_sub: 0.3,
        par: 400.0,
        bw,
        bbr,
        ay
    };


    println!("pp = {:?}", calc_pp(input));

}
