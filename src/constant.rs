//!
//! # Constants
//!
//! This file contains many physical constants useful for scientific applications.
//! They are stored as `f64` to guarantee maximum precision, all are castable to `f32`.
//! 
//! The `constant` modules provides a wide variety of useful mathematical and physical constants. It covers many
//! domains, from quantum mechanics to mathematics, but also thermodynamics, electromagnetism, astrophysics, etc...
//! 
//! Each constant has its uncertainty specified in its description (unless it is defined),
//! and the source is also provided.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Astrophysics and astronomy

/// # $c$ - Speed of light in a vacuum
/// Value is defined
///
/// - Value: $299~792~458$
/// - Unit: $\mathrm{m.s^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?c)
pub const C: f64 = 299_792_458.0;

/// # $G$ - Newtonian gravitational constant
///
/// - Value: $6.674~30\times10^{-11}$
/// - Uncertainty: $0.000~15\times10^{-11}$
/// - Unit: $\mathrm{m^3.kg^{-1}.s^{-2}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?bg)
pub const G: f64 = 6.674_30e-11;

/// # $g$ - Standard earth gravity acceleration
/// Value is defined.
///
/// - Value: $9.806~65$
/// - Unit: $\mathrm{m.s^{-2}}$
/// - Source: [NIST](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication330e2008.pdf)
pub const EARTH_GRAVITY: f64 = 9.806_65;

/// # $M_🜨$ - Earth mass
///
/// - Value: 5.972_167_87×10²⁴ $5.972~167~87\times10^{24}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NASA](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
pub const EARTH_MASS: f64 = 5.972_167_87e24;

/// # $R_🜨$ - Earth radius
///
/// - Value: $6~378~100$
/// - Unit: $\mathrm{m}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/5_English.pdf)
pub const EARTH_RADIUS: f64 = 6_378_100.0;

/// # $G^M_🜨$ - Nominal Earth mass parameter
///
/// - Value: $3.986~004\times10^{14}$
/// - Unit: $\mathrm{m^3.s^{-2}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const EARTH_GM: f64 = 3.986_004e14;

/// # $A_🜨$ - Earth Bond albedo
///
/// - Value: $0.306$
/// - Unit: Dimensionless
pub const EARTH_ALB: f64 = 0.306;

/// # $M_\mathrm{J}$ - Jupiter mass
///
/// - Value: $1.898~124~6\times10^{27}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NASA](https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html)
pub const JUPITER_MASS: f64 = 1.898_124_6e27;

/// # $R_\mathrm{J}$ - Jupiter radius
///
/// - Value: $71~492~000$
/// - Unit: $\mathrm{m}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const JUPITER_RADIUS: f64 = 71_492_000.0;

/// # $G^M_\mathrm{J}$ - Nominal Jupiter mass parameter
///
/// - Value: $1.266~865~3\times10^{17}$
/// - Unit: $\mathrm{m^3.s^{-2}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const JUPITER_GM: f64 = 1.266_865_3e17;

/// # $M_\odot$ - Sun mass
///
/// - Value: $1.988~409~87\times10^{30}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NASA](https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
pub const SUN_MASS: f64 = 1.988_409_87e30;

/// # $R_\odot$ - Sun radius
///
/// - Value: $695 658 000$
/// - Unit: $\mathrm{m}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const SUN_RADIUS: f64 = 695_658_000.0;

/// # $G^M_\odot$ - Nominal Sun mass parameter
///
/// - Value: $1.327~124~4\times10^{20}$
/// - Unit: $\mathrm{m^3.s^{-2}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const SUN_GM: f64 = 1.327_124_4e20;

/// # $T_\odot$ - Sun effective temperature
///
/// - Value: $5772$
/// - Unit: $\mathrm{K}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const SUN_TEFF: f64 = 5772.0;

/// # $f_\odot$ - Sun irradiance
/// 
/// - Value: $1361$
/// - Unit: $\mathrm{W.m^{-2} = kg.s^{-3}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const SUN_IRR: f64 = 1361.0;

/// # $L_\odot$ - Sun luminosity
///
/// - Value: $3.828\times10^{26}$
/// - Unit: $\mathrm{W = kg.m^{2}.s^{-3}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const SUN_L: f64 = 3.828e26;

/// # $L_\circ$ - Absolute bolometric magnitude
///
/// - Value: $3.012~8\times10^{28}$
/// - Unit: $\mathrm{W = kg.m^{2}.s^{-3}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const L_0: f64 = 3.012_75e28;

/// # $M_\mathrm{shift}$ - Absolute bolometric magnitude shift
/// Value is computed based on the apparent bolometric luminosity $L_\circ$
/// and the magnitude definition, leading to:
/// $$
/// M_\mathrm{shift} = \frac{2.5\ln(L_\circ)}{\ln(2) + \ln(5)}
/// $$
/// 
/// - Value: $71.197~425~756~681~473~979~136~000~883~404~700~135~818~216~505~689~428~018~913...$
/// - Unit: Dimensionless
/// - Source: Computed based on [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const ABS_MAG_SHIFT: f64 = 71.197_425_756_681_473_979_136_000_883_404_700_135_818_216_505_689_428_018_913;

/// # $f_\circ$ - Apparent bolometric luminosity
/// 
/// - Value: $2.518\times10^{-8}$
/// - Unit: $\mathrm{W.m^{-2} = kg.s^{-3}}$
/// - Source: [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const F_0: f64 = 2.518_021_002e-8;

/// # $m_\mathrm{shift}$ - Apparent bolometric magnitude shift
/// Value is computed based on the apparent bolometric luminosity $f_\circ$
/// and the magnitude definition, leading to:
/// $$
/// m_\mathrm{shift} = \frac{2.5\ln(f_\circ)}{\ln(2) + \ln(5)}
/// $$
/// 
/// - Value: $-18.997~351~629~757~571~863~584~819~458~253~886~318~389~141~306~029~739~058...$
/// - Unit: Dimensionless
/// - Source: Computed based on [IAU](https://www.iau.org/static/resolutions/IAU2015_English.pdf)
pub const APP_MAG_SHIFT: f64 = -18.997_351_629_757_571_863_584_819_458_253_886_318_389_141_306_029_739_058_11;

/// # $\mathrm{AU}$ - Astronomical unit
/// Value is defined
///
/// - Value: $1.495~978~707\times10^{11}$
/// - Unit: $\mathrm{m}$
/// - Source: [NASA](https://ssd.jpl.nasa.gov/astro_par.html)
pub const AU: f64 = 1.495_978_707e11;

/// # $\mathrm{ly}$ - Light year
/// Value is defined as!
/// $$
/// 1\mathrm{ly} = c\times 1\mathrm{yr}
/// $$
/// Where $c$ is the speed of light and $1\mathrm{yr}$ is a year.
///
/// - Value: $9.460~730~472~580~8\times10^{15}$
/// - Unit: $\mathrm{m}$
/// - Source: Computed from definition
pub const LY: f64 = 9.460_730_472_580_8e15;

/// # $\mathrm{pc}$ - Parsec
/// Value is defined as:
/// $$
/// 1\mathrm{pc} = \frac{1}{\tan(1^{\prime\prime})}\mathrm{AU}
/// $$
/// Where $AU$ is an astronomical unit, and $1^{\prime\prime}$ is an arc-second.
///
/// - Value: $3.085~677~581~491~367~278~913~937~957~796~472\times10^{16}$
/// - Unit: $\mathrm{m}$
/// - Source: Computed from definition
pub const PC: f64 = 3.085_677_581_491_367_278_913_937_957_796_472e16;

/// # $\mathrm{erg}$ - Erg
/// 
/// - Value: $10^{-7}$
/// - Unit: $\mathrm{J}$
/// - Source: [IAU](https://www.iau.org/publications/proceedings_rules/units/)
pub const ERG: f64 = 1e-7;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Universal constants

/// # $h$ - Planck constant
/// Value is defined.
/// 
/// - Value: $6.626~070~15\times10^{-34}$
/// - Unit: $\mathrm{J.s = kg.m^2.s^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?h)
pub const H: f64 = 6.626_070_15e-34;

/// # $\hbar$ - Reduced Planck constant
/// $$
/// \hbar = \frac{h}{2\pi}
/// $$
/// Value is defined, is computed as precise as $2\pi$ is in Rust.
/// 
/// - Value: $1.054~571~817\times10^{-34}$
/// - Unit: $\mathrm{J.s = kg.m^2.s^{-1}}$
/// - Source: Computed from definition [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?hbar)
pub const H_BAR: f64 = H / std::f64::consts::TAU;

/// # $l_\mathrm{P}$ - Planck length
///
/// - Value: $1.616~255\times10^{-35}$
/// - Uncertainty: $0.000~018\times10^{-35}$
/// - Unit: $\mathrm{m}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?plkl)
pub const PLANCK_LENGTH: f64 = 1.616_255e-35;

/// # $m_\mathrm{P}$ - Planck mass
///
/// - Value: $2.176~434\times10^{-8}$
/// - Uncertainty: $0.000~024\times10^{-8}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?plkm)
pub const PLANCK_MASS: f64 = 2.176_434e-8;

/// # $T_\mathrm{P}$ - Planck temperature
///
/// - Value: $1.416~784\times10^{32}$
/// - Uncertainty: $0.000~016\times10^{32}$
/// - Unit: $\mathrm{K}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?plktmp)
pub const PLANCK_TEMP: f64 = 1.416_784e32;

/// # $t_\mathrm{P}$ - Planck time
///
/// - Value: $5.391~247\times10^{-44}$
/// - Uncertainty: $0.000~060\times10^{-44}$
/// - Unit: $\mathrm{s}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?plkt)
pub const PLANCK_TIME: f64 = 5.391_247e-44;

/// # $\epsilon_0$ - Vacuum electric permittivity
///
/// - Value: $8.854~187~812~8\times10^{-12}$
/// - Uncertainty: $0.000~000~001~3\times10^{-12}$
/// - Unit: $\mathrm{F.m^{-1} = s^4.A^2.m^{-3}.kg^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?ep0)
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;

/// # $\mu_0$ - Vacuum magnetic permeability
///
/// - Value: $1.256~637~062~12\times10^{-6}$
/// - Uncertainty: $0.000~000~000~19\times10^{-6}$
/// - Unit: $\mathrm{N.A^{-2} = kg.m.s^{-2}.A^{-2}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?mu0)
pub const MU_0: f64 = 1.256_637_062_12e-6;

/// # $Z_0$ - Characteristic impedance of vacuum
///
/// - Value: $376.730~313~668$
/// - Uncertainty: $0.000~000~057$
/// - Unit: $\mathrm{\Omega = kg.m^2.s^{-3}.A^{-2}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?z0)
pub const Z_0: f64 = 376.730_313_668;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # $\gamma$ - Euler-Mascheroni constant
/// Not to be confused with Euler's number 2.71... (in standard consts of Rust).
/// 
/// - Value: $0.577~215~664~901~532~860~606~512~090~082~402~431~042~159~335~939~92...$
/// - Unit: Dimensionless
pub const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_860_606_512_090_082_402_431_042_159_335_939_92;

/// # $N_\mathrm{A}$ - Avogadro constant
/// Value is defined.
///
/// - Value: $6.022~140~76\times10^{23}$
/// - Unit: $\mathrm{mol^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?na)
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// # $k_\mathrm{B}$ - Boltzmann constant
///
/// - Value: $1.380~649\times10^{-23}$
/// - Unit: $\mathrm{J.K^{-1} = kg.m^2.s^{-2}.K^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?k)
pub const K_B: f64 = 1.380_649e-23;

/// # $e$ - Elementary charge
/// Value is defined.
///
/// - Value: $1.602~176~634\times10^{-19}$
/// - Unit: $\mathrm{C = A.s}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?e)
pub const E: f64 = 1.602_176_634e-19;

/// # $\mathrm{atm}$ - Atmosphere pressure
///
/// - Value: $101~325$
/// - Unit: $\mathrm{Pa = kg.m^{-1}.s^{-2}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?stdatm)
pub const ATM: f64 = 101_325.0;

/// # $\mathrm{bar}$ - Bar pressure unit
/// Value is defined.
/// 
/// - Value: $10^5$
/// - Unit: $\mathrm{Pa = kg.m^{-1}.s^{-2}}$
/// - Source: [IAU](https://www.iau.org/publications/proceedings_rules/units/)
pub const BAR: f64 = 1e5;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Quantum

/// # $a_0$ - Bohr radius
///
/// - Value: $5.291~772~109~03\times10^{-11}$
/// - Uncertainty: $0.000~000~000~80\times10^{-11}$
/// - Unit: $\mathrm{m}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0)
pub const A_0: f64 = 5.291_772_109_03e-11;

/// # $\mu_\mathrm{B}$ - Bohr magneton
///
/// - Value: $9.274~010~078~3\times10^{-24}$
/// - Unit: $\mathrm{J.T^{-1} = kg^2.m^2.s^{-4}.A^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?mub)
pub const BOHR_MAG: f64 = 9.274_010_078_3e-24;

/// # $b$ - Wien wavelength displacement constant
/// Value is defined.
///
/// - Value: $2.897~771~955\times10^{-3}$
/// - Unit: $\mathrm{m.K}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?bwien)
pub const WIEN_B: f64 = 2.897_771_955e-3;

/// # $b'$ - Wien frequency displacement constant
/// Value is defined.
/// 
/// - Value: $5.878~925~757\times10^{10}$
/// - Unit: $\mathrm{Hz.K^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?bpwien)
pub const WIEN_B_FREQ: f64 = 5.878_925_757e10;

/// # $r_e$ - Classical electron radius
///
/// - Value: $2.817~940~326~2×10^{-15}$
/// - Uncertainty: $0.000~000~001~3\times10^{-15}$
/// - Unit: $\mathrm{m}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?re)
pub const E_RADIUS_C: f64 = 2.817_940_326_2e-15;

/// # $\lambda_\mathrm{C}$ - Compton wavelength
///
/// - Value: $2.426~310~238~67\times10^{-12}$
/// - Uncertainty: $0.000~000~000~73\times10^{-12}$
/// - Unit: $\mathrm{m}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?ecomwl)
pub const LAMBDA_COMPTON: f64 = 2.426_310_238_67e-12;

/// # $m_e$ - Electron mass
///
/// - Value: $9.109~383~701~5\times10^{-31}$
/// - Uncertainty: $0.000~000~002~8\times10^{-31}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?me)
pub const ELECTRON_MASS: f64 = 9.109_383_701_5e-31;

/// # $m_\mathrm{n}$ - Neutron mass
///
/// - Value: $1.674~927~498~04\times10^{-27}$
/// - Uncertainty: $0.000~000~000~95\times10^{-27}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?mn)
pub const NEUTRON_MASS: f64 = 1.674_927_498_04e-27;

/// # $m_\mathrm{p}$ - Proton mass
///
/// - Value: $1.672~621~923~69\times10^{-27}$
/// - Uncertainty: $0.000~000~000~51\times10^{-27}$
/// - Unit: $\mathrm{kg}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?mp)
pub const PROTON_MASS: f64 = 1.672_621_923_69e-27;

/// # $R_\infty$ - Rydberg constant
///
/// - Value: $10~973~731.568~160$
/// - Uncertainty: $0.000~021$
/// - Unit: $\mathrm{m^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?ryd)
pub const RYD: f64 = 10_973_731.568_160;

/// # $R$ - Molar gas constant
/// Value is defined.
///
/// - Value: $8.314~462~618$
/// - Unit: $\mathrm{J.mol^{-1}.K^{-1} = kg.m^2.s^{-2}.mol^{-1}.K^{-1}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?r)
pub const R: f64 = 8.314_462_618;

/// # $\sigma$ - Stefan-Boltzmann constant
/// Value is defined.
///
/// - Value: $5.670~374~419\times10^{-8}$
/// - Unit: $\mathrm{W.m^{-2}.K^{-4} = kg.s^{-3}.K^{-4}}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?sigma)
pub const SIGMA_SB: f64 = 5.670_374_419e-8;

/// # $\sigma_\mathrm{e}$ - Thomson scattering cross section
///
/// - Value: $6.652~458~732~1\times10^{-29}$
/// - Unit: $\mathrm{m^2}$
/// - Source: [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?sigmae)
pub const SIGMA_E: f64 = 6.652_458_732_1e-29;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Scales

// Multiples

/// # Yotta prefix
/// 
/// - Value: $10^{24}$
pub const YOTTA: f64 = 1e24;

/// # Zetta prefix
/// 
/// - Value: $10^{21}$
pub const ZETTA: f64 = 1e21;

/// # Exa prefix
/// 
/// Value: $10^{28}$
pub const EXA: f64 = 1e18;

/// # Peta prefix
/// 
/// - Value: $10^{15}$
pub const PETA: f64 = 1e15;

/// # Tera prefix
/// 
/// - Value: $10^{12}$
pub const TERA: f64 = 1e12;

/// # Giga prefix
/// 
/// - Value: $10^{9}$
pub const GIGA: f64 = 1e9;

/// # Mega prefix
/// 
/// - Value: $10^{6}$
pub const MEGA: f64 = 1e6;

/// # Kilo prefix
/// 
/// - Value: $10^{3}$
pub const KILO: f64 = 1e3;

/// # Hecto prefix
///
/// - Value: $10^{2}$
pub const HECTO: f64 = 1e2;

/// # Deca prefix
/// 
/// - Value: $10^{1}$
pub const DECA: f64 = 1e1;

// Sub-multiples

/// # Deci prefix
/// 
/// - Value: $10^{-1}$
pub const DECI: f64 = 1e-1;

/// # Centi prefix
/// 
/// - Value: $10^{-2}$
pub const CENTI: f64 = 1e-2;

/// # Milli prefix
/// 
/// -Value: $10^{-3}$
pub const MILLI: f64 = 1e-3;

/// # Micro prefix
/// 
/// - Value: $10^{-6}$
pub const MICRO: f64 = 1e-6;

/// # Nano prefix
/// 
/// - Value: $10^{-9}$
pub const NANO: f64 = 1e-9;

/// # Pico prefix
/// 
/// - Value: $10^{-12}$
pub const PICO: f64 = 1e-12;

/// # Femto prefix
/// 
/// - Value: $10^{-15}$
pub const FEMTO: f64 = 1e-15;

/// # Atto prefix
/// 
/// - Value: $10^{-18}$
pub const ATTO: f64 = 1e-18;

/// # Zepto prefix
/// 
/// - Value: $10^{-21}$
pub const ZEPTO: f64 = 1e-21;

/// # Yocto prefix
/// 
/// - Value: $10^{-24}$
pub const YOCTO: f64 = 1e-24;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
