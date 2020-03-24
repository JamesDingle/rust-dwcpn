//use crate::modules::time::gen_time_array;

use crate::dwcpn::dwcpn::{InputParams, calc_pp};

mod dwcpn;

fn main() {

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
    };

    println!("pp = {:?}", calc_pp(input));

}
