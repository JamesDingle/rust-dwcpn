use crate::dwcpn::modules::config::NUM_WAVELENGTHS;

const NANO: [f64; NUM_WAVELENGTHS] = [
    0.016     ,  0.0176    ,  0.02426667,  0.0265625 ,  0.02776562,
    0.02896875,  0.03017188,  0.031375  ,  0.03257812,  0.03241667,
    0.03020833,  0.028     ,  0.03135   ,  0.0347    ,  0.03242917,
    0.03015833,  0.0278875 ,  0.02561667,  0.0232619 ,  0.02057143,
    0.01788095,  0.01519048,  0.0125    ,  0.0114    ,  0.0103    ,
    0.0093    ,  0.0083    ,  0.007125  ,  0.00595   ,  0.004775  ,
    0.0036    ,  0.0024    ,  0.0017    ,  0.0007    ,  0.0005    ,
    0.001325  ,  0.00215   ,  0.002975  ,  0.0038    ,  0.00390345,
    0.0040069 ,  0.00411034,  0.00421379,  0.00431724,  0.00443333,
    0.0046    ,  0.0057375 ,  0.006875  ,  0.0080125 ,  0.00915   ,
    0.0102875 ,  0.011425  ,  0.0125625 ,  0.0137    ,  0.0199    ,
    0.02020769,  0.02051538,  0.01826471,  0.01217647,  0.00608824,  0.
];

const PICO: [f64; NUM_WAVELENGTHS] = [
    0.09572727,  0.1053    ,  0.11313333,  0.1197625 ,  0.12609063,
    0.13241875,  0.13874687,  0.145075  ,  0.15140312,  0.15406667,
    0.15123333,  0.1484    ,  0.1401    ,  0.1318    ,  0.12582083,
    0.11984167,  0.1138625 ,  0.10788333,  0.1011    ,  0.0911    ,
    0.0811    ,  0.0711    ,  0.0611    ,  0.0514    ,  0.0417    ,
    0.0351    ,  0.0285    ,  0.025625  ,  0.02275   ,  0.019875  ,
    0.017     ,  0.0155    ,  0.0136    ,  0.0131    ,  0.0122    ,
    0.011525  ,  0.01085   ,  0.010175  ,  0.0095    ,  0.00982759,
    0.01015517,  0.01048276,  0.01081034,  0.01113793,  0.01166667,
    0.013     ,  0.014925  ,  0.01685   ,  0.018775  ,  0.0207    ,
    0.022625  ,  0.02455   ,  0.026475  ,  0.0284    ,  0.0348    ,
    0.03156923,  0.02833846,  0.02329412,  0.01552941,  0.00776471,  0.
];

const MICRO: [f64; NUM_WAVELENGTHS] = [
    0.01518182,  0.0167    ,  0.01686667,  0.0176    ,  0.018475  ,
    0.01935   ,  0.020225  ,  0.0211    ,  0.021975  ,  0.02213333,
    0.02121667,  0.0203    ,  0.01925   ,  0.0182    ,  0.0173875 ,
    0.016575  ,  0.0157625 ,  0.01495   ,  0.0141381 ,  0.01332857,
    0.01251905,  0.01170952,  0.0109    ,  0.0101    ,  0.0093    ,
    0.0087    ,  0.0081    ,  0.0075    ,  0.0069    ,  0.0063    ,
    0.0057    ,  0.0049    ,  0.0043    ,  0.0038    ,  0.0037    ,
    0.0038    ,  0.0039    ,  0.004     ,  0.0041    ,  0.00428966,
    0.00447931,  0.00466897,  0.00485862,  0.00504828,  0.00526667,
    0.0056    ,  0.0066125 ,  0.007625  ,  0.0086375 ,  0.00965   ,
    0.0106625 ,  0.011675  ,  0.0126875 ,  0.0137    ,  0.0165    ,
    0.01480769,  0.01311538,  0.01067647,  0.00711765,  0.00355882,  0.
];

pub fn calc_ac(chl: f64) -> ([f64; NUM_WAVELENGTHS], f64) {

    let mut chlorophyll_absorption: [f64; NUM_WAVELENGTHS] = [0.0; NUM_WAVELENGTHS];
    let mut absorption_sum: f64 = 0.0;

    // coefficients
    // taken from Brewin et al 2011 and 2015
    let cm_pn: f64 = 0.77;
    let s_pn: f64 = 0.94 / cm_pn;
    let cm_p: f64 = 0.13;
    let s_p: f64 = 0.80 / cm_p;

    // compute fractions
    for i in 0..NUM_WAVELENGTHS {
        let mut pico_absorption: f64 = cm_p * ( 1.0 - (-s_p * chl).exp() );
        let mut nano_absorption: f64 = cm_pn * ( 1.0 - ( -s_pn * chl).exp() ) - pico_absorption;
        let mut micro_absorption: f64 = chl - (cm_pn * ( 1.0 - (-s_pn * chl).exp() ));

        // guarantee that no negative absorption gets through as it is phsyically impossible
        if pico_absorption < 0.0 { pico_absorption = 0.0;};
        if nano_absorption < 0.0 { nano_absorption = 0.0;};
        if micro_absorption < 0.0 { micro_absorption = 0.0;};

        chlorophyll_absorption[i] = (PICO[i] * pico_absorption) + (NANO[i] * nano_absorption) + (MICRO[i] * micro_absorption);
        absorption_sum = absorption_sum + chlorophyll_absorption[i];
    }

    let chlorophyll_absorption = chlorophyll_absorption; // de-mutify
    let absorption_mean = absorption_sum / NUM_WAVELENGTHS as f64;

    return (chlorophyll_absorption, absorption_mean);
}