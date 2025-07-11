import numpy as np

def calculate_chi(kh):
    """
    Calculates the stable amplitude factor chi for a given non-dimensional parameter kh.

    The formula is derived from the Stokes flow solution for a viscous layer
    with a harmonic basal shear stress and a rigid, no-slip top boundary.

    chi = - (2*sinh(kh) + kh*cosh(kh)) / (cosh(kh)^2 - (kh)^2)

    Args:
        kh (float): The non-dimensional parameter, product of wavenumber k and layer thickness h.

    Returns:
        float: The calculated amplitude factor chi.
    """
    if kh == 0:
        print("Error: kh cannot be zero as it would lead to division by zero.")
        return None

    # --- Calculate each term in the equation ---
    
    # Value of the input parameter
    val_kh = kh
    
    # Hyperbolic sine of kh
    sinh_kh = np.sinh(kh)
    
    # Hyperbolic cosine of kh
    cosh_kh = np.cosh(kh)
    
    # Numerator of the equation for chi
    # Equation: -(2*sinh(kh) + kh*cosh(kh))
    numerator_val = -(2 * sinh_kh + kh * cosh_kh)
    
    # Denominator of the equation for chi
    # Equation: cosh(kh)^2 - kh^2
    denominator_val = cosh_kh**2 - kh**2
    
    if denominator_val == 0:
        print(f"Error: Denominator is zero for kh = {kh}. Cannot calculate chi.")
        return None

    # Final amplitude factor chi
    chi = numerator_val / denominator_val

    # --- Print the results step-by-step ---
    print(f"Calculating the amplitude factor chi for kh = {val_kh:.4f}")
    print("-" * 50)
    print(f"1. Intermediate values:")
    print(f"   sinh(kh) = {sinh_kh:.4f}")
    print(f"   cosh(kh) = {cosh_kh:.4f}")
    print("-" * 50)
    print(f"2. Final equation components:")
    print(f"   Numerator = -(2 * {sinh_kh:.4f} + {val_kh:.4f} * {cosh_kh:.4f}) = {numerator_val:.4f}")
    print(f"   Denominator = {cosh_kh:.4f}^2 - {val_kh:.4f}^2 = {denominator_val:.4f}")
    print("-" * 50)
    print(f"3. Final Result:")
    print(f"   Amplitude Factor chi = {numerator_val:.4f} / {denominator_val:.4f} = {chi:.4f}")
    
    return chi

# --- User Execution ---
# You can change the value of kh to see how the amplitude factor changes.
# kh represents the ratio of the plate thickness to the wavelength of the stress.
# A value of kh=1 means the thickness is 1/(2*pi) of the wavelength.
kh_value = 1.0
calculate_chi(kh_value)

# Example with a different value
# kh_value_2 = 2.0
# print("\n--- Another Example ---")
# calculate_chi(kh_value_2)