use crate::dwcpn::dwcpn::{calc_pp, InputParams};
use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};

use netcdf;
use pbr::ProgressBar;

mod dwcpn;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // pre calculate bw/bbr/ay arrays for use in the pp_profile calculation later
    // this only needs to be done once for all pixels so I have it here, before we start
    // looping over pixels
    let bw = calculate_bw();
    let bbr = calculate_bbr();
    let ay = calculate_ay();

    // let filename = "/home/jad/Downloads/pp_gku_params_2003_07.nc";
    // let filename = "pp_processing_20100501_30.000_50.000_-140.000_-120.000.nc";
    let filename = "/home/jad/work/pp/istar_testing/atlantic/pp_processing_20150501_49.000_76.000_-47.000_4.000_rust.nc";
    // let filename = "/home/jad/work/pp/czcs_global/pp_processing_19980101_-90.000_90.000_-180.000_180.000.nc";
    // let filename = "/home/jad/work/pp/istar_testing/pp_processing_20100525_-10.000_10.000_120.000_140.000_rust.nc";
    // let filename = "/home/jad/work/pp/new_zenith/global/pp_processing_20100501_-90.000_90.000_-180.000_180.000.nc";

    println!("Processing file: {}", filename);

    let mut ncfile = netcdf::append(&filename)?;

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

    let lat_data = lat.values::<f64>(None, None)?;
    let lon_data = lon.values::<f64>(None, None)?;
    let chl_data = chl.values::<f64>(None, None)?;
    let par_data = par.values::<f64>(None, None)?;
    let bathymetry_data = bathymetry.values::<f64>(None, None)?;
    let pi_alpha_data = pi_alpha.values::<f64>(None, None)?;
    let pi_pmb_data = pi_pmb.values::<f64>(None, None)?;
    let zm_data = zm.values::<f64>(None, None)?;
    let rho_data = rho.values::<f64>(None, None)?;
    let sigma_data = sigma.values::<f64>(None, None)?;

    let mut pp_data = vec![9969209968386869000000000000000000000.0; lat.len() * lon.len()];
    let mut spectral_i_star_mean_data = vec![9969209968386869000000000000000000000.0; lat.len() * lon.len()];
    let mut euphotic_depth_data =
        vec![9969209968386869000000000000000000000.0; lat.len() * lon.len()];

    println!("{:?}", lat.len());
    println!("{:?}", lon.len());

    let count = lat.len() * lon.len();
    let mut pb = ProgressBar::new(count as u64);

    for y in 0..lat.len() {
        for x in 0..lon.len() {
            pb.inc();

            if (chl_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (par_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (bathymetry_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (pi_alpha_data[[y, x]] == -999.0)
                || (pi_pmb_data[[y, x]] == -999.0)
                || (zm_data[[y, x]] == -999.0)
                || (rho_data[[y, x]] == -999.0)
                || (sigma_data[[y, x]] == -999.0)
            {
                continue;
            }

            let input = InputParams {
                lat: lat_data[y],
                lon: lon_data[x],
                z_bottom: bathymetry_data[[y, x]],
                iday: 135,
                // iday: 165,
                alpha_b: pi_alpha_data[[y, x]],
                pmb: pi_pmb_data[[y, x]],
                z_m: zm_data[[y, x]],
                chl: chl_data[[y, x]],
                rho: rho_data[[y, x]],
                sigma: sigma_data[[y, x]],
                cloud: 0.0,
                yel_sub: 0.3,
                par: par_data[[y, x]],
                bw,
                bbr,
                ay,
            };

            let result = calc_pp(input);
            if result.0 < 10000.0 {
                pp_data[y * lon.len() + x] = result.0;
                euphotic_depth_data[y * lon.len() + x] = result.1;
                spectral_i_star_mean_data[y * lon.len() + x] = result.2;
            }
        } // x loop
    } // y loop

    // let mut ncfile2 = netcdf::append(&filename)?;
    let mut pp_var = ncfile.variable_mut("pp").unwrap();
    &pp_var.put_values(&pp_data, None, None);

    let mut euph_var = ncfile.variable_mut("euphotic_depth").unwrap();
    &euph_var.put_values(&euphotic_depth_data, None, None);

    let mut spectral_i_star_var = ncfile.variable_mut("spectral_i_star_mean").unwrap();
    &spectral_i_star_var.put_values(&spectral_i_star_mean_data, None, None);
    Ok(())
}
