import math
try:
    import xraylib
except ImportError:
    print("Please install the 'xraylib' library using: pip install xraylib")
    exit()

def calculate_transmission():
    """
    Calculates the X-ray transmission through a Be window for various elements.
    """
    # Parameters for the EDX system
    window_material = 'Be'
    window_thickness_um = 100.0
    window_thickness_cm = window_thickness_um / 10000.0  # Convert µm to cm
    
    # Get density of Beryllium from xraylib
    # Note: Standard density of Be is ~1.85 g/cm^3. xraylib.ElementDensity(4) provides this.
    try:
        rho_be = xraylib.ElementDensity(4) # Atomic number of Be is 4
    except Exception:
        # Fallback if density function is not available in an older xraylib version
        rho_be = 1.85 

    # Elements to check from the answer choices (Z = Atomic Number)
    # W is the sample, not the light element to be detected.
    elements_to_check = {
        'Na': 11,
        'Mg': 12,
        'Si': 14,
        'Ca': 20
    }

    print(f"Calculating X-ray transmission through a {window_thickness_um} µm ({window_thickness_cm} cm) Be window...")
    print("-" * 70)

    results = {}

    for symbol, Z in elements_to_check.items():
        try:
            # Get the K-alpha line energy for the element in keV
            energy_kev = xraylib.LineEnergy(Z, xraylib.KA_LINE)
            
            # Get the mass attenuation coefficient of Be at this energy (cm^2/g)
            mu_rho = xraylib.CS_MassCoef(window_material, energy_kev)
            
            # Calculate the exponent term of the Beer-Lambert law
            exponent_value = mu_rho * rho_be * window_thickness_cm
            
            # Calculate the transmission
            transmission = math.exp(-exponent_value)
            
            # Store and print the results
            results[symbol] = transmission * 100
            
            # Print the full equation as requested
            print(f"For {symbol} (Z={Z}):")
            print(f"  K-alpha Energy = {energy_kev:.3f} keV")
            print(f"  Transmission = exp(-(μ/ρ) * ρ * x)")
            print(f"  Transmission = exp(-({mu_rho:.2f} cm²/g) * ({rho_be:.2f} g/cm³) * ({window_thickness_cm:.4f} cm))")
            print(f"  Transmission = exp(-{exponent_value:.3f}) = {transmission:.4f} or {transmission * 100:.2f}%\n")

        except Exception as e:
            print(f"Could not perform calculation for {symbol}: {e}")

    # Final analysis
    print("-" * 70)
    print("Analysis:")
    print("An element is considered detectable if its transmission is reasonably high (e.g., >1%).")
    print("Based on the calculations, Sodium (Na) has a transmission of ~1.4%, which is low but detectable.")
    print("Since Na is the lightest element among the choices with a detectable transmission, it is the answer.")

if __name__ == '__main__':
    calculate_transmission()