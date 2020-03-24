use crate::dwcpn::modules::config::NUM_WAVELENGTHS;

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

// correct factors for diffuse irradiance for 5 wavelengths at 7 zenith angles
const CORRECTION: [[f64; 5]; 7] = [
    [1.11, 1.04, 1.15, 1.12, 1.32],
    [1.13, 1.05, 1.00, 0.96, 1.12],
    [1.18, 1.09, 1.00, 0.96, 1.07],
    [1.24, 1.11, 0.99, 0.94, 1.02],
    [1.46, 1.24, 1.06, 0.99, 1.10],
    [1.70, 1.34, 1.07, 0.96, 0.90],
    [2.61, 1.72, 1.22, 1.04, 0.80]
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

    for wl in 0..TRANSMITTANCE_WL_COUNT {
        let wld = TRANSMITTANCE_WAVELENGTHS[wl] / 1000.0; // NOTE: This is repeated a few times, worth a single pre-compute?
        if wl < 10 {
            ta[wl] = (-BETA1 * wld.powf(-ALPHA1) * airmass).exp();
        } else {
            ta[wl] = (-BETA2 * wld.powf(-ALPHA2) * airmass).exp();
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


// given the calculated zenith angle (DEGREES), look up the index of the correction
// factors we need to use
fn find_lut_index(zen_d: f64) -> usize {
    for i in 0..CORRECTION_ZEN_LOOKUP.len() {
        if zen_d < CORRECTION_ZEN_LOOKUP[i] {
            return i
        }
    }
    return CORRECTION_ZEN_LOOKUP.len() - 1
}

// take the correction factors for 5 wavelengths from the LUT and interpolate to the
// number of wavelengths we are using in the transmittance calculations (24 at time of writing)
fn interpolate_correction_factor(zen_d: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {

    let lut_index: usize = find_lut_index(zen_d);

    let c: [f64; 5] = CORRECTION[lut_index];
    let cm1: [f64; 5] = CORRECTION[lut_index - 1];

    let fraction: f64 = (zen_d - CORRECTION_ZEN_LOOKUP[lut_index - 1]) / (CORRECTION_ZEN_LOOKUP[lut_index] - CORRECTION_ZEN_LOOKUP[lut_index - 1]);

    let mut c1: [f64; 5] = [0.0; 5];
    for i in 0..c.len() {
        c1[i] = (c[i] - cm1[i]) * fraction + cm1[i];
    }
    let c1 = c1;

    let mut c2: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    let mut cinc: f64 = 0.0;
    let mut l: usize = 0;

    c2[0] = c1[0];

    for l1 in 1..4 {
        cinc = (c1[l1 - 1] - c1[l1]) / 5.0;
        for l2 in 0..5 {
            l += 1;
            c2[l] = c2[l - 1] - cinc;
        }
    }

    // for the remaining wavelengths 550 -> 710nm
    let cdif: f64 = c1[4] - c1[3];

    let wldif = TRANSMITTANCE_WAVELENGTHS[16] - TRANSMITTANCE_WAVELENGTHS[TRANSMITTANCE_WL_COUNT - 1];

    println!("{:?}", wldif);

    [0.0; TRANSMITTANCE_WL_COUNT]

}

fn compute_diffuse_irradiance(
    zen_r: f64, zen_d: f64,
    direct_irradiance: [f64; NUM_WAVELENGTHS],
    ro_s: [f64; NUM_WAVELENGTHS],
    t_aerosol: [f64; NUM_WAVELENGTHS],
    t_ozone: [f64; NUM_WAVELENGTHS],
    t_rayleigh: [f64; NUM_WAVELENGTHS],
    tu: f64,
    t_water_vapour: [f64; NUM_WAVELENGTHS]
) -> [f64; NUM_WAVELENGTHS] {

    let correction_matrix: [f64; TRANSMITTANCE_WL_COUNT] = interpolate_correction_factor(zen_d);

    [0.0; NUM_WAVELENGTHS]
}