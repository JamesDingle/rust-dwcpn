//use crate::modules::time::gen_time_array;

use crate::dwcpn::dwcpn::{InputParams, calc_pp};
use crate::dwcpn::modules::pp_profile::{calculate_bw, calculate_bbr, calculate_ay};

use netcdf;
use pbr::ProgressBar;

mod dwcpn;


// You can also modify a Variable inside an existing netCDF file
// open it in read/write mode
// let mut file = netcdf::append("crabs2.nc")?;
// // get a mutable binding of the variable "crab_coolness_level"
// let mut var = file.variable_mut("crab_coolness_level").unwrap();
//
// let data : Vec<i32> = vec![100; 10];
// // write 5 first elements of the vector `data` into `var` starting at index 2;
// var.put_values(&data[..5], Some(&[2]), Some(&[5]));
// // Change the first value of `var` into '999'
// var.put_value(999.0f32, Some(&[0]));

fn main() {

    // pre calculate bw/bbr/ay arrays for use in the pp_profile calculation later
    // this only needs to be done once for all pixels so I have it here, before we start
    // looping over pixels
    let bw = calculate_bw();
    let bbr = calculate_bbr();
    let ay = calculate_ay();

    // let filename = "/home/jad/Downloads/pp_gku_params_2003_07.nc";
    // let filename = "pp_processing_20100501_30.000_50.000_-140.000_-120.000.nc";
    let filename = "/tmp/ppws/pp_processing_19980101_-50.000_-30.000_160.000_180.000.nc";

    let ncfile = netcdf::append(&filename).unwrap();

    let lat = &ncfile.variable("lat").unwrap();
    let lon = &ncfile.variable("lon").unwrap();
    let chl = &ncfile.variable("chlor_a").unwrap();
    let par = &ncfile.variable("par").unwrap();
    let bathymetry = &ncfile.variable("bathymetry").unwrap();
    let pi_alpha = &ncfile.variable("PI_alpha").unwrap();
    let pi_pmb = &ncfile.variable("PI_pmb").unwrap();
    let zm = &ncfile.variable("zm").unwrap();
    let rho = &ncfile.variable("rho").unwrap();
    let sigma = &ncfile.variable("sigma").unwrap();

    // let mut pp_var: [f64; 241*241] = [0.0; 241 * 241];

    println!("{:?}", lat.len());
    println!("{:?}", lon.len());
    let mut ncfile2 = netcdf::append(&filename).unwrap();
    let mut pp:netcdf::variable::VariableMut = ncfile2.variable_mut("pp").unwrap();

    let count = lat.len() * lon.len();
    let mut pb = ProgressBar::new(count as u64);

    for y in 0..lat.len() {
        for x in 0..lon.len() {
            pb.inc();
            let input = InputParams{
                lat: lat.value(Some(&[y])).unwrap(),
                lon: lon.value(Some(&[x])).unwrap(),
                z_bottom: bathymetry.value(Some(&[y,x])).unwrap(),
                iday: 100,
                alpha_b: pi_alpha.value(Some(&[y,x])).unwrap(),
                pmb: pi_pmb.value(Some(&[y,x])).unwrap(),
                z_m: zm.value(Some(&[y,x])).unwrap(),
                chl: chl.value(Some(&[y,x])).unwrap(),
                rho: rho.value(Some(&[y,x])).unwrap(),
                sigma: sigma.value(Some(&[y,x])).unwrap(),
                cloud: 0.0,
                yel_sub: 0.3,
                par: par.value(Some(&[y,x])).unwrap(),
                bw,
                bbr,
                ay
            };

            if (input.chl == 9969209968386869000000000000000000000.0) ||
                (input.par == 9969209968386869000000000000000000000.0) ||
                (input.z_bottom == 9969209968386869000000000000000000000.0) ||
                (input.alpha_b == 9969209968386869000000000000000000000.0) ||
                (input.pmb == 9969209968386869000000000000000000000.0) ||
                (input.z_m == 9969209968386869000000000000000000000.0) ||
                (input.rho == 9969209968386869000000000000000000000.0) ||
                (input.sigma == 9969209968386869000000000000000000000.0) {

                pp.put_value(9969209968386869000000000000000000000.0, Some(&[y,x])).unwrap();
                continue
            }

            let result = calc_pp(input);

            if result.0 < 4000.0 {
                // println!("pp = {:?}", result);
                // pp_var[x * y + y] = result.0;
                // pp.put_value(result.0, Some(&[y,x])).unwrap();
                pp.put_value(result.0, Some(&[y,x])).unwrap();

            } else {
                pp.put_value(9969209968386869000000000000000000000.0, Some(&[y,x])).unwrap();
                // println!("test");
            }

        } // x loop

    } // y loop

    // pb.finish_print("done");

    // let mut ncfile = ncfile.borrow_mut();

    // let mut pp_var = ncfile.variable_mut("pp").unwrap();



}

