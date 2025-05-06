"""
# BHFallback.jl

This file contains functions for modeling black hole fallback mass ejection. These functions calculate the remaining mass of a black hole after a fraction of its mass is ejected, with optional restrictions based on the initial mass.

## Overview
- `arbitraryEjectaBH`: Calculates the remaining mass of a black hole after an arbitrary fraction of its mass is ejected.
- `restrictedEjectaBH`: Calculates the remaining mass of a black hole with restrictions on the ejected fraction based on the initial mass.

These functions are useful for astrophysical simulations involving black hole mass loss and fallback scenarios.
"""


"""
    arbitraryEjectaBH(m2_i, frac)

    This function calculates the remaining mass of a black hole after a fraction of its mass is ejected. The fraction is specified by the `frac` parameter.

    # Arguments:
    - `m2_i`: Initial mass of the black hole.
    - `frac`: Fraction of the mass to be ejected.

    # Output:
    - Returns the remaining mass of the black hole after the ejection.
"""
function arbitraryEjectaBH(m2_i, frac)
    return m2_i*(1-frac)
end

"""
TODO: Fix this
    restrictedEjectaBH(m2_i, frac; Ma=10, Mb=15, max_frac=1.0, min_frac=0.1)
    This function calculates the remaining mass of a black hole after a fraction of its mass is ejected, with restrictions based on the initial mass. The fraction is specified by the `frac` parameter.
    The function applies different limits on the fraction of mass that can be ejected based on the initial mass of the black hole. If the initial mass is between `Ma` and `Mb`, the maximum fraction is adjusted. If the initial mass is greater than `Mb`, the minimum fraction is applied.
    # Arguments:
    - `m2_i`: Initial mass of the black hole.
    - `frac`: Fraction of the mass to be ejected.
    - `Ma`: Lower mass limit for the adjustment (default is 10).
    - `Mb`: Upper mass limit for the adjustment (default is 15).
    - `max_frac`: Maximum fraction of mass that can be ejected (default is 1.0).
    - `min_frac`: Minimum fraction of mass that can be ejected (default is 0.1).
    # Output:
    - Returns the remaining mass of the black hole after the ejection, adjusted based on the initial mass.
"""
function restrictedEjectaBH(m2_i, frac; Ma=10, Mb=15, max_frac=1.0, min_frac=0.1)
    frac_limit = max_frac
    if m2_i > Ma && m2_i < Mb
        frac_limit = max_frac + (min_frac-max_frac)*(m2_i-Ma)/(Mb-Ma)
    elseif m2_i > Mb
        frac_limit = min_frac
    end
    return m2_i*(1-frac*frac_limit)
end
