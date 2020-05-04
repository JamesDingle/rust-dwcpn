pub const TIMESTEPS: usize = 12;

pub const DEPTH_PROFILE_COUNT: usize = 500;
pub const DEPTH_PROFILE_STEP: f64 = 0.5;
pub const DEPTH_PROFILE_START: f64 = 0.0;

pub const NUM_WAVELENGTHS: usize = 61;

pub const WAVELENGTHS: [f64; NUM_WAVELENGTHS] = [
    400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0,
    465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0,
    530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0,
    595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0,
    660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0
];

pub const AW: [f64; NUM_WAVELENGTHS] = [
    0.00663, 0.00530, 0.00473, 0.00444,
    0.00454, 0.00478, 0.00495, 0.00530,
    0.00635, 0.00751, 0.00922, 0.00962,
    0.00979, 0.01011, 0.0106, 0.0114,
    0.0127, 0.0136, 0.0150, 0.0173,
    0.0204, 0.0256, 0.0325, 0.0396,
    0.0409, 0.0417, 0.0434, 0.0452,
    0.0474, 0.0511, 0.0565, 0.0596,
    0.0619, 0.0642, 0.0695, 0.0772,
    0.0896, 0.1100, 0.1351, 0.1672,
    0.2224, 0.2577, 0.2644, 0.2678,
    0.2755, 0.2834, 0.2916, 0.3012,
    0.3108, 0.325, 0.340, 0.371,
    0.410, 0.429, 0.439, 0.448,
    0.465, 0.486, 0.516, 0.559, 0.624
];