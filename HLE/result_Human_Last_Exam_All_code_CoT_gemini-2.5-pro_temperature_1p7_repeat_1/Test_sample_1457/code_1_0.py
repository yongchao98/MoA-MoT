import math
import xraylib

def solve_edx_detection_limit():
    """
    Calculates the lightest element detectable in an EDX system with a thick Be window.

    This script determines the detectability of elements by calculating the transmission
    of their characteristic K-alpha X-rays through a 100-micron Beryllium window.
    An element is detectable if its X-ray transmission is significant (e.g., > 1%).
    """

    # --- Setup Parameters ---
    # Beryllium (Be) window properties
    be_atomic_number = 4
    be_density_g_cm3 = 1.85
    be_thickness_cm = 100 * 1e-4  # 100 µm = 0.01 cm

    # Candidate elements from the multiple-choice question
    # Format: {Symbol: Atomic Number (Z)}
    # We are looking for the lightest element, so we sort them by Z.
    candidate_elements = {
        'Na': 11,
        'Mg': 12,
        'Si': 14,
        'Ca': 20
    }

    print("--- EDX Light Element Detection Analysis ---")
    print(f"Detector Window: {be_thickness_cm * 1e4} µm Beryllium (Density: {be_density_g_cm3} g/cm^3)\n")

    results = []

    # --- Calculations ---
    for symbol, Z in sorted(candidate_elements.items(), key=lambda item: item[1]):
        try:
            # Get K-alpha X-ray energy in keV for the element
            # xraylib.LineEnergy returns eV, so we convert to keV
            e_kev = xraylib.LineEnergy(Z, xraylib.KA_LINE) / 1000.0

            # Get the mass attenuation coefficient (mu) of Be at this energy
            # The unit is cm^2/g
            mu_be_cm2_g = xraylib.CS_Total(be_atomic_number, e_kev)

            # Calculate the exponent term for Beer-Lambert Law: -(mu * rho * t)
            exponent = -mu_be_cm2_g * be_density_g_cm3 * be_thickness_cm
            
            # Calculate transmission percentage
            transmission = math.exp(exponent) * 100

            results.append({
                'symbol': symbol,
                'Z': Z,
                'energy': e_kev,
                'transmission': transmission,
                'mu': mu_be_cm2_g,
                'exponent': exponent
            })

        except Exception as e:
            print(f"Could not calculate for {symbol}: {e}")

    # --- Output Results ---
    print("Results (sorted by Atomic Number Z):")
    print("-" * 65)
    print(f"{'Element':<10} {'Z':<5} {'Kα Energy (keV)':<18} {'Transmission (%)':<20}")
    print("-" * 65)

    lightest_detectable = None

    for res in results:
        print(f"{res['symbol']:<10} {res['Z']:<5} {res['energy']:<18.3f} {res['transmission']:<20.2f}")
        # A common threshold for detectability is >1% transmission.
        # We are looking for the LIGHTEST element that meets this criterion.
        if res['transmission'] > 1.0 and lightest_detectable is None:
            lightest_detectable = res
    
    print("-" * 65)
    
    if lightest_detectable:
        print("\n--- Conclusion ---")
        print(f"The lightest element on the list is Sodium (Na, Z=11).")
        print("To see if it's detectable, we check its Kα X-ray transmission through the window.")
        print(f"Calculation: Transmission = exp(-μ * ρ * t)")
        print(f"For Sodium (Na):")
        print(f"  Kα Energy (E) = {lightest_detectable['energy']:.3f} keV")
        print(f"  Be Mass Attenuation (μ) at E = {lightest_detectable['mu']:.2f} cm²/g")
        print(f"  Be Density (ρ) = {be_density_g_cm3} g/cm³")
        print(f"  Be Thickness (t) = {be_thickness_cm} cm")
        print(f"  Exponent = -({lightest_detectable['mu']:.2f} * {be_density_g_cm3} * {be_thickness_cm}) = {lightest_detectable['exponent']:.3f}")
        print(f"  Transmission = exp({lightest_detectable['exponent']:.3f}) * 100% = {lightest_detectable['transmission']:.2f}%")
        print(f"\nA transmission of {lightest_detectable['transmission']:.2f}% is significant enough to be detected.")
        print(f"Therefore, Sodium (Na) is the lightest element among the choices that can be seen.")
    else:
        print("None of the elements have significant transmission through the window.")


if __name__ == '__main__':
    # Ensure xraylib is installed, or provide instructions
    try:
        import xraylib
    except ImportError:
        print("Error: The 'xraylib' library is required. Please install it using 'pip install xraylib'")
    else:
        solve_edx_detection_limit()

<<<B>>>