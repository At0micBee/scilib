//!
//! # Methods for astronomy
//!

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{     // Using std lib constants
    PI,                     // Pi
    FRAC_PI_2               // Pi / 2
};

use std::fmt::{             // Formatter display
    Display,                // The display itself
    Result as DRes          // The associated result
};

use super::constant;
use super::coordinate::spherical::Spherical;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Radec coordinate system
///
/// ## Definition
/// Right ascension and declination of the object in the sky. The values are stored as `f64` internally, and in radians.
/// The distance is kept in meters internally, for homogeneousness.
/// 
/// ## Example
/// ```
/// # use scilib::astronomy::Radec;
/// let radec = Radec { ra: 0.3, dec: 1.2, dist_earth: None };
/// assert_eq!(radec.ra, 0.3);
/// assert_eq!(radec.dec, 1.2);
/// assert_eq!(radec.dist_earth, None);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Radec {
    /// Right ascension
    pub ra: f64,
    /// Declination
    pub dec: f64,
    /// Distance from earth
    pub dist_earth: Option<f64>
}

/// # Display
///
/// Shows both ra and dec in degrees, which are more easily readable.
impl Display for Radec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> DRes {
        write!(f, "ra={}° :: dec={}°", self.ra.to_degrees(), self.dec.to_degrees())?;

        // If there is a distance, we print it as well
        match self.dist_earth {
            Some(d) => write!(f, " :: distance={}ly", d / constant::LY)?,
            None => {}
        }
        Ok(())
    }
}

/// Implementing required methods
impl Radec {

    /// # Radec coordinates from radians
    /// Generates the struct from angles given in radians.
    /// 
    /// ## Inputs
    /// - `ra`: right ascension, in radians
    /// - `dec`:  declination, in radians
    /// 
    /// Returns a new Radec struct.
    /// 
    /// ## Example
    /// ```
    /// # use std::f64::consts::{ FRAC_PI_3, FRAC_PI_6 };
    /// # use scilib::astronomy::Radec;
    /// let coord1 = Radec::from_rad(FRAC_PI_6, FRAC_PI_3);
    /// let coord2 = Radec { ra: FRAC_PI_6, dec: FRAC_PI_3, dist_earth: None };
    /// assert!((coord1.ra - coord2.ra).abs() < 1.0e-15);
    /// assert!((coord1.dec - coord2.dec).abs() < 1.0e-15);
    /// ```
    pub fn from_rad<T, U>(ra: T, dec: U) -> Self
    where T: Into<f64>, U: Into<f64> {
        Self {
            ra: ra.into(),
            dec: dec.into(),
            dist_earth: None
        }
    }

    /// # Radec coordinates from angles in degrees
    /// Generates the struct from angles given in degrees.
    ///
    /// ## Inputs
    /// - `ra`: right ascension, in degrees
    /// - `dec`:  declination, in degrees
    /// 
    /// Returns a new Radec struct.
    /// 
    /// ## Example
    /// ```
    /// # use std::f64::consts::{ FRAC_PI_3, FRAC_PI_6 };
    /// # use scilib::astronomy::Radec;
    /// let coord1 = Radec::from_degree(30.0, 60.0);
    /// let coord2 = Radec { ra: FRAC_PI_6, dec: FRAC_PI_3, dist_earth: None };
    /// assert!((coord1.ra - coord2.ra).abs() < 1.0e-15);
    /// assert!((coord1.dec - coord2.dec).abs() < 1.0e-15);
    /// ```
    pub fn from_degree<T, U>(ra: T, dec: U) -> Self
    where T: Into<f64>, U: Into<f64> {
        Self {
            ra: ra.into().to_radians(),
            dec: dec.into().to_radians(),
            dist_earth: None
        }
    }

    /// # Specifying known distance
    /// Appends the distance to the Radec struct.
    /// 
    /// ## Input
    /// - `distance`: distance of the object from earth, in meters ($m$)
    /// 
    /// Returns the initial struct with the distance appended.
    /// 
    /// ## Example
    /// ```
    /// # use scilib::constant;
    /// # use scilib::astronomy::Radec;
    /// let mut coord = Radec::from_degree(30.0, 12.3);
    /// coord.specify_distance(4.403 * constant::LY);
    /// assert_eq!(coord.ra, 30.0_f64.to_radians());
    /// assert_eq!(coord.dec, 12.3_f64.to_radians());
    /// assert_eq!(coord.dist_earth, Some(4.403 * constant::LY));
    /// ```
    pub fn specify_distance<T>(&mut self, distance: T)
    where T: Into<f64> {
        self.dist_earth = Some(distance.into());
    }

    /// # Angular separation computation
    ///
    /// ## Definition
    /// We use the Vincenty formula to compute the angular separation, as it should remain stable no matter
    /// the separations and positions of the objects. The resulting separation will be in the range [0, $\pi$];
    ///
    /// ## Inputs
    /// - `self`: the first object in Radec
    /// - `other`: the second object in Radec
    /// 
    /// Returns the angular separation in radians ($rad$).
    /// 
    /// ## Example
    /// ```
    /// # use scilib::astronomy::Radec;
    /// let c1 = Radec::from_rad(1.0, -1.2);
    /// let c2 = Radec::from_rad(-2.02, 0.13);
    /// assert!((c1.separation(&c2) - 2.0685709648870154).abs() < 1.0e-15);
    /// ```
    pub fn separation(&self, other: &Self) -> f64 {

        // Difference in right ascension
        let da: f64 = self.ra - other.ra;

        //Computing the two top terms
        let t1: f64 = (other.dec.cos() * da.sin()).powi(2);
        let t2: f64 = (self.dec.cos() * other.dec.sin() - self.dec.sin() * other.dec.cos() * da.cos()).powi(2);

        // Computing the bottom term
        let b: f64 = self.dec.sin() * other.dec.sin() + self.dec.cos() * other.dec.cos() * da.cos();

        let res: f64 = ((t1 + t2).sqrt() / b).atan();

        if res < 0.0 {
            res + PI
        } else {
            res
        }
    }

    /// # Distance in space between two objects
    /// 
    /// ## Definition
    /// Computes the actual distance in space between two objects, using the spherical coordinates
    /// 
    /// ## Inputs
    /// - `self`: the first object in Radec
    /// - `other`: the second object in Radec
    /// 
    /// Returns the distance between the two objects.
    /// 
    /// ## Example
    /// ```
    /// # use std::f64::consts::SQRT_2;
    /// # use scilib::astronomy::Radec;
    /// let mut c1 = Radec::from_degree(0, 0);
    /// c1.specify_distance(1.0);
    /// 
    /// let mut c2 = Radec::from_degree(0, 90);
    /// c2.specify_distance(1.0);
    /// 
    /// assert!((c1.distance_with(c2) - SQRT_2).abs() < 1.0e-12);
    /// ```
    pub fn distance_with(self, other: Self) -> f64 {

        let obj_1: Spherical = self.into();
        let obj_2: Spherical = other.into();

        obj_1.distance(obj_2)
    }
}

/// # Conversion to spherical coordinates
impl Into<Spherical> for Radec {
    fn into(self) -> Spherical {
        Spherical {
            r: self.dist_earth.unwrap_or_default(),
            theta: self.ra,
            phi: FRAC_PI_2 - self.dec
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Apparent magnitude
/// 
/// ## Definition
/// Uses the standard reference luminosity to compute the apparent magnitude of an object.
/// $$
/// m = -2.5\log\left( \frac{F}{F_\circ} \right)
/// $$
/// Where $F$ is the irradiance of the object, computed with $F = I(L, d)$
/// and $F_\circ$ is the apparent bolometric luminosity.
/// 
/// ## Inputs
/// - `lum`: luminosity of the object ($L$), in watts ($W$)
/// - `dist`: distance from the observer ($d$), in meters ($m$)
/// 
/// Returns the apparent magnitude of the object, dimensionless.
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::apparent_mag;
/// let alpha_b = 0.4 * constant::SUN_L;        // Alpha centauri B, roughly
/// let distance = 4.403 * constant::LY;        // Seen from Earth
/// let res = apparent_mag(alpha_b, distance);
/// assert!((res - 1.36).abs() <= 0.03);
/// ```
pub fn apparent_mag(lum: f64, dist: f64) -> f64 {
    // We avoid a division by using the pre-computed value
    -2.5 * irradiance(lum, dist).log10() + constant::APP_MAG_SHIFT
}

/// # Absolute magnitude
/// 
/// ## Definition
/// Same approach that the apparent magnitude, but set at a distance of 10pc.
/// $$
/// M = -2.5\log\left( \frac{L}{L_\circ} \right)
/// $$
/// Where $L_\circ$ is the absolute bolometric magnitude.
/// 
/// ## Inputs
/// - `lum`: luminosity of the object ($L$), in watts ($W$)
/// 
/// Returns the magnitude of the object, dimensionless.
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::absolute_mag;
/// let abs = absolute_mag(constant::SUN_L);    // Sun absolute magnitude
/// assert!((abs - 4.74).abs() <= 0.02);        // As defined IAU
/// ```
pub fn absolute_mag(lum: f64) -> f64 {
    // We avoid a division by using the pre-computed value
    -2.5 * lum.log10() + constant::ABS_MAG_SHIFT
}

/// # Distance modulus
/// 
/// ## Definition
/// Computes the distance of an object based on its apparent and absolute magnitudes:
/// $$
/// d = 10^{1 + \frac{m-M}{5}}
/// $$
/// Where $m$ is the apparent magnitude and $M$ is the absolute magnitude.
/// 
/// ## Inputs
/// - `m_app`: apparent magnitude ($m$), dimensionless
/// - `m_abs`: absolute magnitude ($M$), dimensionless
/// 
/// Returns the distance in parsecs ($pc$).
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::{ apparent_mag, absolute_mag, distance_mod };
/// let app = apparent_mag(constant::SUN_L, 10.0 * constant::PC);
/// let abs = absolute_mag(constant::SUN_L);
/// let res = distance_mod(app, abs) / constant::PC;
/// assert!((res - 10.0).abs() <= 1e-5);
/// ```
pub fn distance_mod(m_app: f64, m_abs: f64) -> f64 {
    10.0_f64.powf(1.0 + (m_app - m_abs) / 5.0) * constant::PC
}

/// # Equilibrium temperature
/// 
/// ## Definition
/// Knowing the stellar effective temperature, its radius, the distance from the star
/// and the albedo of the considered planet, we can compute the equilibrium temperature
/// found at the planet, computed with:
/// $$
/// T_\mathrm{eq} = T_\mathrm{star} \sqrt{\frac{R_\mathrm{star}}{2d}} (1 - A)
/// $$
/// Where $T_\mathrm{star}$ and $R_\mathrm{star}$ are the effective temperature of the host star and its radius,
/// $d$ is the distance from the star and $A$ is the albedo of the planet.
/// 
/// ## Inputs
/// - `star_t`: star effective temperature ($T_\mathrm{star}$), in kelvins ($K$)
/// - `star_rad`: star radius ($R_\mathrm{star}$), in meters ($m$)
/// - `dist`: distance to host star ($d$), in meters ($m$)
/// - `albedo`: albedo of the planet ($A$), dimensionless
/// 
/// Returns the temperature in kelvins ($K$).
pub fn equilibrium_temperature(star_t: f64, star_rad: f64, dist: f64, albedo: f64) -> f64 {
    // Following the formula
    star_t * (star_rad / (2.0 * dist)).sqrt() * (1.0 - albedo).sqrt()
}

/// # Object irradiance
/// 
/// ## Definition
/// Computes the irradiance at a certain distance from an object, using its luminosity.
/// The `luminosity` is in Watts and the `distance` is in meters, result is in `W.m-2`.
/// The formula used is:
/// $$
/// I(L, d) = \frac{L}{4\pi d^2}
/// $$
/// Where $L$ is the luminosity and $d$ is the distance.
///
/// ## Inputs
/// - `luminosity`: luminosity of the object ($L$), in watts ($W$)
/// - `distance`: distance form the object ($d$), in meters ($m$)
/// 
/// Returns the irradiance in watt per square meters ($W.m^{-2}$)
/// 
/// ## Example
/// ```rust
/// # use scilib::constant;
/// # use scilib::astronomy;
/// // Computing the solar irradiance at the Earth
/// let p_earth: f64 = astronomy::irradiance(constant::SUN_L, constant::AU);
/// assert!((p_earth - 1361.0).abs() < 1.0);
///
/// // Computing the irradiance from the sun at Alpha-Centauri
/// let p: f64 = astronomy::irradiance(constant::SUN_L, 4.3 * constant::LY);
/// assert!((p - 1.84e-8).abs() < 0.1e-8);
/// ```
pub fn irradiance(luminosity: f64, distance: f64) -> f64 {
    // Following the formula
    luminosity / (4.0 * PI * distance.powi(2))
}

/// # Energy received by an object
/// 
/// ## Definition
/// Computes the received energy by an object of a given surface, at a known distance.
/// Makes use of the `irradiance` function to compute the surfacic power output by the object,
/// it is defined as:
/// $$
/// E = I(L, d) S
/// $$
/// Where $I(L, d)$ is the irradiance function and $S$ is the surface of the object.
/// 
/// ## Inputs
/// - `luminosity`: luminosity of the object ($L$), in watts ($W$)
/// - `distance`: distance form the object ($d$), in meters ($m$)
/// - `surface`: exposed surface of the object ($S$), in square meters ($m^2$)
/// 
/// Returns the irradiance of the object in watts ($W$).
pub fn received_energy(luminosity: f64, distance: f64, surface: f64) -> f64 {
    // Following the formula
    irradiance(luminosity, distance) * surface
}

/// # Planetary luminosity
/// 
/// ## Definition
/// Computes the luminosity of a planet, based on the received irradiance and the albedo:
/// $$
/// L_\mathrm{P} = E_\mathrm{received} A = I(L, d) S A
/// $$
/// Where $S$ is the surface, $A$ the albedo, and $E_\mathrm{received} = I(L, d) S$ is the energy
/// received by the planet.
/// 
/// ## Inputs
/// - `albedo`: albedo of the object ($A$), dimensionless
/// - `received`: power received by the object ($I(L, d)$), in watts ($W$)
/// 
/// Returns the luminosity in watts ($W$).
pub fn planet_luminosity(albedo: f64, received: f64) -> f64 {
    received * albedo
}

/// # Luminosity using Stefan-Boltzmann
/// 
/// ## Definition
/// Computes the expected luminosity of a star using the Stefan-Boltzmann constant, from:
/// $$
/// L = 4\pi\sigma R^2 T^4
/// $$
/// 
/// ## Inputs
/// - `radius`: radius of the star ($R$), in meters ($m$)
/// - `temperature`: effective temperature of the star ($T$), in kelvins ($K$)
/// 
/// Returns the luminosity of the star in watts ($W$).
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::{ luminosity, absolute_mag };
/// let sun = luminosity(constant::SUN_RADIUS, constant:: SUN_TEFF);
/// assert!((sun - constant::SUN_L).abs() / constant::SUN_L <= 1.0e-3);
/// 
/// let computed = absolute_mag(sun);
/// let theory = absolute_mag(constant::SUN_L);
/// assert!((computed - theory).abs() / theory <= 1.0e-4);
/// ```
pub fn luminosity(radius: f64, temperature: f64) -> f64 {
    4.0 * PI * constant::SIGMA_SB * radius.powi(2) * temperature.powi(4)
}

/// # Rayleigh criterion
/// ## Definition
/// The Rayleigh criterion qualifies the minimum separation of two points that can be resolved
/// as two distinct points by an instrument. It equates simply as:
/// $$
/// \epsilon = 1.22\frac{\lambda}{D}
/// $$
/// Where $\lambda$ is the considered wavelength and $D$ is the diameter of the instrument.
/// 
/// ## Inputs
/// - `lambda`: considered wavelength ($\lambda$), in meters ($m$)
/// - `d`: the diameter of the instrument ($d$), in meters ($m$)
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::rayleigh_criterion;
/// let wave: f64 = 500.0 * constant::NANO;     // Greenish light
/// let diam: f64 = 8.2;                        // VLT diameter
/// let criterion = rayleigh_criterion(wave, diam);
/// assert!((criterion - 7.43902439024e-8).abs() < 1.0e-8)
/// ```
pub fn rayleigh_criterion(lambda: f64, d: f64) -> f64 {
    1.22 * lambda / d
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Orbital speed of an object
/// ## Definition
/// Computation for the velocity of an object on an elliptical orbit. The definition is:
/// $$
/// v = \sqrt{GM\left( \frac{2}{R} - \frac{1}{a} \right)}
/// $$
/// Where $G$ is the gravitational constant, $M$ is the mass of the orbited body, $R$ is the current distance
/// of the object from the center of mass of the orbited body and $a$ is the semi-major axis of the orbit.
/// 
/// ## Inputs
/// - `mass`: the mass of the center body ($M$), in kilograms ($kg$)
/// - `r`: distance from the center of focal point ($R$), in meters ($m$)
/// - `a`: semi-major axis of the orbit ($a$), in meters ($m$)
/// 
/// Returns $v$, the velocity of an object at a given altitude.
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::orbital_speed;
/// let rad = constant::EARTH_RADIUS + 421.0e3;     // ISS perigee
/// let semi = 6795.0e3;                            // ISS semi-major axis
/// let v = orbital_speed(constant::EARTH_MASS, rad, semi);
/// assert!((v - 7654.416466312335).abs() < 1.0e-10);
/// ```
pub fn orbital_speed(mass: f64, r: f64, a: f64) -> f64 {
    ( constant::G * mass * (2.0 / r - 1.0 / a) ).sqrt()
}

/// # Impact parameter $b$
/// ## Definition
/// The impact parameter for an exoplanet orbiting a star is defined as:
/// $$
/// b = \frac{a\cdot\cos(i)}{R_s}
/// $$
/// 
/// ## Inputs
/// - `a`: semi-major axis ($a$), in meters ($m$)
/// - `radius_star`: the radius of the host star ($R_s$), in meters ($m$)
/// - `i`: inclination of the planet's orbit ($i$), in degrees
/// 
/// Return $b$, the impact parameter of the planet.
/// 
/// ## Example
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::impact_parameter;
/// let r_star: f64 = 0.834 * constant::SUN_RADIUS;
/// let a: f64 = 0.07697 * constant::AU;
/// let i: f64 = 88.7;
/// let b: f64 = impact_parameter(a, r_star, i);
/// 
/// assert!((b - 0.474).abs() <= 0.025);
/// ```
pub fn impact_parameter(a: f64, radius_star: f64, i: f64) -> f64 {
    i.to_radians().cos() * a / radius_star
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
