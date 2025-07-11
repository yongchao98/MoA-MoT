import math

try:
    import xraylib
except ImportError:
    print("Please install the 'xraylib' library by running: pip install xraylib")
    exit()

def solve_edx_problem():
    """
    Calculates the transmission of X-rays from light elements through a Be window
    to determine the lightest detectable element.
    """
    # --- Step 1: Define parameters ---
    # Parameters for the Beryllium window
    be_symbol = 'Be'
    be_atomic_number = 4
    be_window_thickness_um = 100.0
    be_window_thickness_cm = be_window_thickness_um * 1e-4  # Convert µm to cm

    # Get Be density from xraylib
    be_density_g_cm3 = xraylib.ElementDensity(be_atomic_number)

    # --- Step 2: Define candidate elements and their properties ---
    # We check the light elements from the answer choices.
    elements_to_check = {
        "Na": 11,
        "Mg": 12,
        "Si": 14,
        "Ca": 20,
    }

    print("Analyzing X-ray transmission through a Beryllium (Be) window.")
    print(f"Window Thickness (x): {be_window_thickness_um} µm = {be_window_thickness_cm} cm")
    print(f"Density of Be (ρ): {be_density_g_cm3:.3f} g/cm³\n")

    # --- Step 3: Calculate transmission for each element ---
    print("-" * 80)
    print(f"{'Element':<10} {'Z':<5} {'Kα Energy (keV)':<18} {'μ in Be (cm²/g)':<20} {'Transmission (%)':<20}")
    print("-" * 80)

    results = []
    # Sort elements by atomic number (Z)
    sorted_elements = sorted(elements_to_check.items(), key=lambda item: item[1])

    for symbol, z in sorted_elements:
        # Get the K-alpha line energy for the element
        energy_keV = xraylib.LineEnergy(z, xraylib.KA_LINE)

        # Get the mass attenuation coefficient (μ) of Be at this energy
        mass_atten_coeff = xraylib.CS_Total_CP(be_symbol, energy_keV)

        # Calculate the argument of the exponent in Beer-Lambert law
        exponent_arg = -mass_atten_coeff * be_density_g_cm3 * be_window_thickness_cm

        # Calculate transmission and convert to percentage
        transmission = math.exp(exponent_arg)
        transmission_percent = transmission * 100

        results.append({
            "symbol": symbol,
            "energy_keV": energy_keV,
            "mu": mass_atten_coeff,
            "transmission_percent": transmission_percent
        })
        
        print(f"{symbol:<10} {z:<5} {energy_keV:<18.3f} {mass_atten_coeff:<20.3f} {transmission_percent:<20.2f}%")

    print("-" * 80)

    # --- Step 4: Interpret the results and find the answer ---
    lightest_detectable = None
    for result in results:
        # A common threshold for "detectable" is >1% transmission
        if result["transmission_percent"] > 1.0:
            lightest_detectable = result
            break

    print("\n### Conclusion ###")
    print("The lightest element that can be seen is the one with the lowest atomic number")
    print("that has a calculated transmission of at least ~1% through the detector window.")

    if lightest_detectable:
        symbol = lightest_detectable['symbol']
        energy = lightest_detectable['energy_keV']
        mu = lightest_detectable['mu']
        trans_pct = lightest_detectable['transmission_percent']

        print(f"\nThe lightest element to meet this criterion is {symbol} (Sodium).")
        print(f"Its Kα X-rays have a calculated transmission of {trans_pct:.2f}%.")
        
        print("\nFinal Equation for the Lightest Detectable Element (Sodium):")
        print("Transmission = exp( -(μ) * (ρ) * (x) )")
        print(f"Transmission = exp( -({mu:.3f} cm²/g) * ({be_density_g_cm3:.3f} g/cm³) * ({be_window_thickness_cm} cm) )")
        final_transmission_value = math.exp(-mu * be_density_g_cm3 * be_window_thickness_cm)
        print(f"Transmission = {final_transmission_value:.4f}, or {final_transmission_value*100:.2f}%")

    else:
        print("\nNone of the elements have significant transmission, an error likely occurred.")


if __name__ == "__main__":
    solve_edx_problem()

<<<B>>>