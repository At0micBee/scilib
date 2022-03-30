//!
//! # Constants
//!
//! This file contains many physical constants useful for scientific applications.
//! They are stored as `f64` to guarantee maximum precision, most are castable to `f32`.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Astrophysics and astronomy

/// # Speed of light in a vacuum
///
/// - Value: 299 792 458
/// - Unit: `m.s-1`
/// - Source: [NASA](https://ssd.jpl.nasa.gov/astro_par.html)
pub const C: f64 = 299_792_458.0;

/// # Newtonian gravitational constant
///
/// - Value: 6.674_30×10⁻¹¹
/// - Uncertainty: 0.000 15×10⁻¹¹
/// - Unit: `m3.kg-1.s-2`
/// - Source: [NASA](https://ssd.jpl.nasa.gov/astro_par.html)
pub const G: f64 = 6.674_30e-11;

/// # Standard earth gravity acceleration
///
/// - Value: 9.806 65
/// - Unit: `m.s-2`
pub const EARTH_GRAVITY: f64 = 9.806_65;

/// # Earth mass
///
/// - Value: 5.972_167_87×10²⁴
/// - Unit: `kg`
pub const EARTH_MASS: f64 = 5.972_167_87e24;

/// # Earth radius
///
/// - Value: 6 378 100
/// - Unit: `m`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const EARTH_RADIUS: f64 = 6_378_100.0;

/// # Nominal Earth mass parameter
///
/// - Value: 3.986 004×10¹⁴
/// - Unit: `m3.s-2`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const EARTH_GM: f64 = 3.986_004e14;

/// # Earth Bond albedo
///
/// - Value: 0.306
/// - Unit: /
pub const EARTH_ALB: f64 = 0.306;

/// # Jupiter mass
///
/// - Value: 1.898_124_6×10²⁷
/// - Unit: `kg`
pub const JUPITER_MASS: f64 = 1.898_124_6e27;

/// # Jupiter radius
///
/// - Value: 71 492 000
/// - Unit: `m`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const JUPITER_RADIUS: f64 = 71_492_000.0;

/// # Nominal Jupiter mass parameter
///
/// - Value: 1.266_865_3×10¹⁷
/// - Unit: `m3.s-2`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const JUPITER_GM: f64 = 1.266_865_3e17;

/// # Sun mass
///
/// - Value: 1.988 409 87×10³⁰
/// - Unit: `kg`
pub const SUN_MASS: f64 = 1.988_409_87e30;

/// # Sun radius
///
/// - Value: 695 700 000
/// - Unit: `m`
pub const SUN_RADIUS: f64 = 695_700_000.0;

/// # Nominal Sun mass parameter
///
/// - Value: 1.327 124 4×10²⁰
/// - Unit: `m3.s-2`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const SUN_GM: f64 = 1.327_124_4e20;

/// # Sun effective temperature
///
/// - Value: 5772
/// - Unit: `K`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const SUN_TEFF: f64 = 5772.0;

/// # Sun irradiance
/// 
/// - Value: 1361
/// - Unit: `W.m-2` = `kg.s−3`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const SUN_IRR: f64 = 1361.0;

/// # Sun luminosity
///
/// - Value: 3.828×10²⁶
/// - Unit: `W` = `kg.m2.s−3`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const SUN_L: f64 = 3.828e26;

/// # Absolute bolometric magnitude
///
/// - Value: 3.012_8×10²⁸
/// - Unit: `W` = `kg.m2.s−3`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const L0: f64 = 3.012_8e28;

/// # Apparent bolometric luminosity
/// 
/// - Value: 2.518×10⁻⁸
/// - Unit: `W.m-2` = `kg.s−3`
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
pub const F_0: f64 = 2.518_021_002e-8;

/// # Astronomical unit
///
/// - Value: 1.495 978 71×10¹¹
/// - Unit: `m`
/// - Source: [NASA](https://ssd.jpl.nasa.gov/astro_par.html)
pub const AU: f64 = 1.495_978_70700e11;

/// # Light year
///
/// - Value: 9.460 7×10¹⁵
/// - Unit: `m`
pub const LY: f64 = 9.460_7e15;

/// # Parsec
///
/// - Value: 3.085 677 58×10¹⁶
/// - Unit: `m`
pub const PC: f64 = 3.085_677_58e16;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Universal constants

/// # Planck constant
///
/// - Value: 6.626 070 15×10⁻³⁴
/// - Unit: `J.s` = `kg.m2.s-1`
pub const H: f64 = 6.626_070_15e-34;

/// # Reduced Planck constant (h/2pi)
///
/// - Value: 1.054_571_817×10⁻³⁴
/// - Unit: `J.s` = `kg.m2.s-1`
pub const H_BAR: f64 = 1.054_571_817e-34;

/// # Planck length
///
/// - Value: 1.616 255×10⁻³⁵
/// - Uncertainty: 0.000 018×10⁻³⁵
/// - Unit: `m`
pub const PLANCK_LENGTH: f64 = 1.616_255e-35;

/// # Planck mass
///
/// - Value: 2.176 434×10⁻⁸
/// - Uncertainty: 0.000 024×10⁻⁸
/// - Unit: `kg`
pub const PLANCK_MASS: f64 = 2.176_434e-8;

/// # Planck temperature
///
/// - Value: 1.416_784×10³²
/// - Uncertainty: 0.000 016×10³²
/// - Unit: `K`
pub const PLANCK_TEMP: f64 = 1.416_784e32;

/// # Planck time
///
/// - Value: 5.391 247×10⁻⁴⁴
/// - Uncertainty: 0.000 060×10⁻⁴⁴
/// - Unit: `s`
pub const PLANCK_TIME: f64 = 5.391_247e-44;

/// # Vacuum electric permittivity
///
/// - Value: 8.854 187 812 8×10⁻¹²
/// - Uncertainty: 0.000 000 0013×10⁻¹²
/// - Unit: `F.m-1` = `s4.A2.m-3.kg-1`
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;

/// # Vacuum magnetic permeability
///
/// - Value: 1.256 637 062 12×10⁻⁶
/// - Uncertainty: 0.000 000 000 19×10⁻⁶
/// - Unit: `N.A-2` = `kg.m.s-2.A-2`
pub const MU_0: f64 = 1.256_637_062_12e-6;

/// # Characteristic impedance of vacuum
///
/// - Value: 376.730 313 668
/// - Uncertainty: 0.000 000 057
/// - Unit: `Ohm` = `kg.m2.s-3.A-2`
pub const Z_0: f64 = 376.730_313_668;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defined

/// # Euler-Mascheroni constant (γ)
///
/// - Value: 0.577 215 664 901 532 860 606 512 090 082 402 431 042 159 335 939 92
/// - Unit: /
///
/// Not to be confused with Euler's number 2.71... (in standard consts of Rust).
pub const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_860_606_512_090_082_402_431_042_159_335_939_92;

/// # Avogadro constant
///
/// - Value: 6.022 140 76×10²³
/// - Unit: `mol-1`
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// # Boltzmann constant
///
/// - Value: 1.380_649×10⁻²³
/// - Unit: `J.K-1` = `kg.m2.s-2.k-1`
pub const K_B: f64 = 1.380_649e-23;

/// # Elementary charge
///
/// - Value: 1.602 176 634×10⁻¹⁹
/// - Unit: `C` = `A.s`
pub const E: f64 = 1.602_176_634e-19;

/// # Atmosphere pressure
///
/// - Value: 101 325
/// - Unit: `Pa` = `[kg.m-1.s-2]`
pub const ATM: f64 = 101_325.0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Quantum

/// # Bohr radius
///
/// - Value: 5.291 772 109 03×10⁻¹¹
/// - Uncertainty: 0.000 000 000 80×10⁻¹¹
/// - Unit: `m`
pub const A_0: f64 = 5.291_772_109_03e-11;

/// # Bohr magneton
///
/// - Value: 9.274 010 08×10⁻²⁴
/// - Unit: `J.T-1` = `kg2.m2.s-4.A-1`
pub const BOHR_MAG: f64 = 9.274_010_08e-24;

/// # Wien wavelength
///
/// - Value: 0.002 897 771 96
/// - Unit: `m.K`
pub const WIEN_B: f64 = 0.002_897_771_96;

/// # Classical electron radius
///
/// - Value: 2.817 940 326 2×10⁻¹⁵
/// - Uncertainty: 0.000 000 0013×10⁻¹⁵
/// - Unit: `m`
pub const E_RADIUS_C: f64 = 2.817_940_3262e-15;

/// # Compton wavelength
///
/// - Value: 2.426 310 238 67×10⁻¹²
/// - Uncertainty: 0.000 000 000 73×10⁻¹²
/// - Unit: `m`
pub const LAMBDA_COMPTON: f64 = 2.426_310_238_67e-12;

/// # Electron mass
///
/// - Value: 9.109 383 701 5×10⁻³¹
/// - Uncertainty: 0.000 000 0028×10⁻³¹
/// - Unit: `kg`
pub const ELECTRON_MASS: f64 = 9.109_383_7015e-31;

/// # Neutron mass
///
/// - Value: 1.674 927 498 04×10⁻²⁷
/// - Uncertainty: 0.000 000 000 95×10⁻²⁷
/// - Unit: `kg`
pub const NEUTRON_MASS: f64 = 1.674_927_498_04e-27;

/// # Proton mass
///
/// - Value: 1.672 621 923 69×10⁻²⁷
/// - Uncertainty: 0.000 000 000 51×10⁻²⁷
/// - Unit: `kg`
pub const PROTON_MASS: f64 = 1.672_621_923_69e-27;

/// # Rydberg constant
///
/// - Value: 10 973 731.568 160
/// - Uncertainty: 0.000 021
/// - Unit: `m-1`
pub const RYD: f64 = 10_973_731.568_160;

/// # Gas constant
///
/// - Value: 8.314 462 618
/// - Unit: `J.mol-1.K-1` = `kg.m2.s-2.mol-1.K-1`
pub const R: f64 = 8.314_462_618;

/// # Stefan Boltzmann constant
///
/// - Value: 5.670_374_419×10⁻⁸
/// - Unit: `W.m-2.K-4` = `kg.s−3.K-4`
pub const SIGMA_SB: f64 = 5.670_374_419e-8;

/// # Thomson scattering cross section
///
/// - Value: 6.652 458 73×10⁻²⁹
/// - Unit: `m2`
pub const SIGMA_T: f64 = 6.652_458_73e-29;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Scales

// Multiples

/// # 1e24
/// 1 000 000 000 000 000 000 000 000
pub const YOTTA: f64 = 1e24;

/// # 1e21
/// 1 000 000 000 000 000 000 000
pub const ZETTA: f64 = 1e21;

/// # 1e18
/// 1 000 000 000 000 000 000
pub const EXA: f64 = 1e18;

/// # 1e15
/// 1 000 000 000 000 000
pub const PETA: f64 = 1e15;

/// # 1e12
/// 1 000 000 000 000
pub const TERA: f64 = 1e12;

/// # 1e9
/// 1 000 000 000
pub const GIGA: f64 = 1e9;

/// # 1e6
/// 1 000 000
pub const MEGA: f64 = 1e6;

/// # 1e3
/// 1 000
pub const KILO: f64 = 1e3;

/// # 1e2
/// 100
pub const HECTO: f64 = 1e2;

/// # 1e1
/// 10
pub const DECA: f64 = 1e1;

// Submultiples

/// # 1e-1
/// 0.1
pub const DECI: f64 = 1e-1;

/// # 1e-2
/// 0.01
pub const CENTI: f64 = 1e-2;

/// # 1e-3
/// 0.001
pub const MILLI: f64 = 1e-3;

/// # 1e-6
/// 0.000 001
pub const MICRO: f64 = 1e-6;

/// # 1e-9
/// 0.000 000 001
pub const NANO: f64 = 1e-9;

/// # 1e-12
/// 0.000 000 000 001
pub const PICO: f64 = 1e-12;

/// # 1e-15
/// 0.000 000 000 000 001
pub const FEMTO: f64 = 1e-15;

/// # 1e-18
/// 0.000 000 000 000 000 001
pub const ATTO: f64 = 1e-18;

/// # 1e-21
/// 0.000 000 000 000 000 000 001
pub const ZEPTO: f64 = 1e-21;

/// # 1e-24
/// 0.000 000 000 000 000 000 000 001
pub const YOCTO: f64 = 1e-24;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
