"""
# Constants.jl

This file defines physical and astronomical constants used throughout the SideKicks.jl package.

## Overview
The constants include:
- Gravitational constant (`cgrav`)
- Masses of the Sun and Earth (`m_sun`, `m_earth`)
- Astronomical unit (`au`)
- Radius of the Sun (`r_sun`)
- Time units (`day`, `year`)
- Conversion factors (`km`, `km_per_s`, `degree`)

All constants are defined in CGS (centimeter-gram-second) units unless otherwise noted.
"""

export cgrav, m_sun, m_earth, au, r_sun, km, day, year, km_per_s, degree

# Gravitational constant (cm^3 g^-1 s^-2), CODATA 2018
const cgrav = 6.67430e-8

# GM product for Earth and Sun (cm^3 s^-2), IAU Resolution B3
const mu_sun = 1.3271244e26
const mu_earth = 3.986004e20

# Mass of the Sun and Earth (g)
const m_sun = mu_sun / cgrav
const m_earth = mu_earth / cgrav

# Astronomical unit (cm), IAU 2009 system of astronomical constants
const au = 1.49597870700e13

# Radius of the Sun (cm), Schmutz & Kosovichev (2008)
const r_sun = 6.9566e10

# Conversion factors
const km = 1.0e5          # Kilometers to centimeters
const km_per_s = km       # Kilometers per second (for convenience)

# Time units
const day = 24.0 * 3600.0 # Seconds in a day
const year = 365.25 * day # Seconds in a year

# Angular conversion
const degree = pi / 180   # Degrees to radians