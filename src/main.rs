use crate::dwcpn::dwcpn::{calc_pp, ModelInputs, ModelSettings};
use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};

use clap::{Arg, App, value_t};
use netcdf;
use indicatif::{ProgressBar, ProgressStyle};

mod dwcpn;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // Argument parsing
    let args = App::new("DWCPN Primary Production Model")
        .version("0.1.0")
        .arg(Arg::with_name("inputfile")
            .short("i")
            .long("inputfile")
            .help("location of netcdf file to run the model on")
            .required(true)
            .takes_value(true)
        )
        .arg(Arg::with_name("jday")
            .short("j")
            .long("jday")
            .help("Julian day of year (a.k.a. Ordinal)")
            .required(true)
            .takes_value(true)
        )
        .arg(Arg::with_name("mld")
            .short("m")
            .long("mixed-layer-depth")
            .help("If enabled, limits depth profile to mixed layer depth with surface chlor_a concentration propagated throughout")
            .required(false)
            .takes_value(false)
        )
        .get_matches();

    let filename = args.value_of("filename").unwrap();
    let jday = value_t!(args.value_of("jday"), u16).unwrap();

    let settings = ModelSettings {
        mld_only: args.is_present("mld")
    };

    // pre calculate bw/bbr/ay arrays for use in the pp_profile calculation later
    // this only needs to be done once for all pixels so I have it here, before we start
    // looping over pixels
    let bw = calculate_bw();
    let bbr = calculate_bbr();
    let ay = calculate_ay();

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
    let mld = &ncfile.variable("mld").unwrap();
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
    let mld_data = mld.values::<f64>(None, None)?;
    let rho_data = rho.values::<f64>(None, None)?;
    let sigma_data = sigma.values::<f64>(None, None)?;

    let mut pp_data = vec![9969209968386869000000000000000000000.0; lat.len() * lon.len()];
    let mut spectral_i_star_mean_data = vec![9969209968386869000000000000000000000.0; lat.len() * lon.len()];
    let mut euphotic_depth_data =
        vec![9969209968386869000000000000000000000.0; lat.len() * lon.len()];

    println!("{:?}", lat.len());
    println!("{:?}", lon.len());

    let count = lat.len() * lon.len();
    let pb = ProgressBar::new(count as u64);
    pb.set_draw_delta(100);
    pb.set_style(ProgressStyle::default_bar()
        .template("{msg} {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {percent}% [{pos:>7}/{len:7} @ {per_sec}] (ETA: {eta})")
        .progress_chars("#>-"));

    // for y in 0..lat.len() {
    //     for x in 0..lon.len() {
    for y in 0..lat.len() {
        for x in 0..lon.len() {
            pb.inc(1);

            if (chl_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (par_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (bathymetry_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (mld_data[[y, x]] == 9969209968386869000000000000000000000.0)
                || (pi_alpha_data[[y, x]] == -999.0)
                || (pi_pmb_data[[y, x]] == -999.0)
                || (zm_data[[y, x]] == -999.0)
                || (rho_data[[y, x]] == -999.0)
                || (sigma_data[[y, x]] == -999.0)
            {
                continue;
            }

            let input = ModelInputs {
                lat: lat_data[y],
                lon: lon_data[x],
                z_bottom: bathymetry_data[[y, x]],
                iday: jday,
                // iday: 165,
                alpha_b: pi_alpha_data[[y, x]],
                pmb: pi_pmb_data[[y, x]],
                z_m: zm_data[[y, x]],
                mld: mld_data[[y, x]],
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

            let result = calc_pp(&input, &settings);
            if result.0 < 10000.0 && result.0 > 0.0 {
                pp_data[y * lon.len() + x] = result.0;
                euphotic_depth_data[y * lon.len() + x] = result.1;
                spectral_i_star_mean_data[y * lon.len() + x] = result.2;
            }
        } // x loop
    } // y loop

    pb.finish_with_message("Complete!");

    // let mut ncfile2 = netcdf::append(&filename)?;
    let mut pp_var = ncfile.variable_mut("pp").unwrap();
    &pp_var.put_values(&pp_data, None, None);

    let mut euph_var = ncfile.variable_mut("euphotic_depth").unwrap();
    &euph_var.put_values(&euphotic_depth_data, None, None);


    let mut spectral_i_star_var = ncfile.variable_mut("spectral_i_star_mean").unwrap();
    &spectral_i_star_var.put_values(&spectral_i_star_mean_data, None, None);
    Ok(())
}
