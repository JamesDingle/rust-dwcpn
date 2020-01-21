use crate::dwcpn::modules::config::NUM_WAVELENGTHS;
use test::NamePadding::PadNone;

const TRANSMITTANCE_WL_COUNT: usize = 24;
const TRANSMITTANCE_WAVELENGTHS: [f64; TRANSMITTANCE_WL_COUNT] = [
    400.000, 410.000, 420.000, 430.000, 440.000, 450.000, 460.000, 470.000,
    480.000, 490.000, 500.000, 510.000, 520.000, 530.000, 540.000, 550.000,
    570.000, 593.000, 610.000, 630.000, 656.000, 667.600, 690.000, 710.000
];

// ozone absorption coefficients
const OZONE_ABS: [f64; TRANSMITTANCE_WL_COUNT] = [
    0.000, 0.000, 0.000, 0.000, 0.000, 0.003, 0.006, 0.009,
    0.014, 0.021, 0.030, 0.040, 0.048, 0.063, 0.075, 0.095,
    0.120, 0.119, 0.132, 0.120, 0.065, 0.060, 0.028, 0.018
];


// Water vapour absorption is set to 2.0
const W:f64 = 2.0;
// water vapour absorption coefficients
const WATER_VAPOUR_ABS: [f64; TRANSMITTANCE_WL_COUNT] = [
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
    0.000, 0.075, 0.000, 0.000, 0.000, 0.000, 0.016, 0.0125
];

// Aerosol transmittance coefficients
const ALPHA1: f64 = 1.0274;
const BETA1: f64 = 0.1324;
const ALPHA2: f64 = 1.206;
const BETA2: f64 = 0.117;

// correct factors for diffuse irradiance for 5 wavelengths and 7 zenith angles
const CORRECTION: [[f64; 7]; 5] = [
    [1.11, 1.13, 1.18, 1.24, 1.46, 1.70, 2.61],
    [1.04, 1.05, 1.09, 1.11, 1.24, 1.34, 1.72],
    [1.15, 1.00, 1.00, 0.99, 1.06, 1.07, 1.22],
    [1.12, 0.96, 0.96, 0.94, 0.99, 0.96, 1.04],
    [1.32, 1.12, 1.07, 1.02, 1.10, 0.90, 0.80]
];

const CORRECTION_ZEN_LOOKUP: [f64; 7] = [0., 37., 48.19, 60., 70., 75., 80.];

fn compute_airmass(zen_r: f64, zen_d: f64) -> f64 {
    let mut airmass: f64;

    airmass = 1.0 / (zen_r.cos() + 0.15 * (93.885 - zen_d).powf(-1.253));

    if airmass < 1.0 {
        airmass = 1.0;
    }

    airmass
}

fn compute_rayleigh(airmass: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut tr: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        let wld = TRANSMITTANCE_WAVELENGTHS[w] / 1000.0;
        tr[w] = (-airmass / wld.powf(4.0) * (115.6406 - 1.335 / wld.powf(2.0)));
    }

    tr
}

fn compute_aerosol_transmittance(airmass: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut ta: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        let wld = TRANSMITTANCE_WAVELENGTHS[w] / 1000.0; // NOTE: This is repeated a few times, worth a single pre-compute?
        if i < 10 {
            ta[w] = (-BETA1 * wld.powf(-ALPHA1) * airmass).exp();
        } else {
            ta[w] = (-BETA2 * wld.powf(-ALPHA2) * airmass).exp();
        }
    }
    ta
}

fn compute_water_vapour_transmittance(airmass: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut tw: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        let wld = TRANSMITTANCE_WAVELENGTHS[w] / 1000.0; // NOTE: This is repeated a few times, worth a single pre-compute?
        tw[w] = (-0.3285 * WATER_VAPOUR_ABS[w] * (W + (1.42 - W) / 2.0) * airmass / (1.0 + 20.07 * WATER_VAPOUR_ABS[w] * airmass).powf(0.45)).exp()
    }
    tw
}

fn compute_ozone_transmittance(airmass: f64, zenith_r: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut to: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        let em0: f64 = 35.0 / (1224.0 * (zenith_r.cos()).powf(2.0) + 1.0).powf(0.5);
        to[w] = (-OZONE_ABS[w] * 0.03 * em0).exp();
    }

    to

}

fn compute_tu(airmass: f64) -> f64 {
    (-1.41 * 0.15 * airmass / (1.0 + 118.3 * 0.15 * airmass).powf(0.45)).exp()
}

// double* compute_air_albedo(wavelength_array_t *wl_array, double* ta, double* to, double* tr, double tu, double* tw){
//    double* ro_s = (double*)malloc(sizeof(double) * wl_array->count);
//
//    int i;
//    for (i = 0; i < wl_array->count; ++i){
//        ro_s[i] = to[i] * tw[i] * (ta[i] * (1.0 - tr[i]) * 0.5 + tr[i] * (1 - ta[i]) * 0.22 * 0.928);
//    }
//    ro_s[23] = ro_s[23] * tu;
//    return ro_s;
//}

fn compute_air_albedo(
    aerosol: [f64; TRANSMITTANCE_WL_COUNT],
    ozone: [f64; TRANSMITTANCE_WL_COUNT],
    rayleigh: [f64; TRANSMITTANCE_WL_COUNT],
    water_vapour: [f64; TRANSMITTANCE_WL_COUNT],
    tu: f64
) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut rho_s: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0.. TRANSMITTANCE_WL_COUNT {
        rho_s[w] = ozone[w] * water_vapour[w] * (aerosol[w] * (1.0 - rayleigh[w]) * 0.5 + rayleigh[w] * (1.0 - aerosol[w]) * 0.22 * 0.928);
    }

    rho_s[23] = rho_s[23] * tu;
    rho_s
}

//double* compute_diffuse_irradiance(wavelength_array_t *wl_array, double zenith_d, double zenith_r, const double* direct, const double* ro_s, const double* ta, const double* to, const double* tr, double tu, const double* tw){

//fn compute_diffuse_irradiance(
//    zen_r: f64, zen_d: f64,
//    direct_irradiance: [f64; NUM_WAVELENGTHS],
//    ro
//) -> [f64; NUM_WAVELENGTHS] {
//
//}