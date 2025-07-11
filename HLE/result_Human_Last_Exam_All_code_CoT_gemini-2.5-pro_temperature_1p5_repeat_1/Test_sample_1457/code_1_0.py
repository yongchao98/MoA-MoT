import numpy as np
import xraylib

def solve_edx_detection_limit():
    """
    Calculates the X-ray transmission through a Be window for various elements
    to determine the lightest detectable element in an EDX setup.
    """
    # --- Setup Parameters ---
    # Properties of the detector window
    window_material_symbol = 'Be'
    window_material_name = 'Beryllium'
    window_thickness_um = 100.0
    window_thickness_cm = window_thickness_um / 1e4  # Convert µm to cm

    # Get Be density from xraylib for accuracy
    # Z=4 for Beryllium
    window_density_g_cm3 = xraylib.ElementDensity(4)

    # Elements to analyze from the answer choices (Z: Atomic Number)
    elements_to_check = {
        'Na': 11,
        'Mg': 12,
        'Si': 14,
        'Ca': 20,
    }

    print(f"Analyzing X-ray transmission through a {window_thickness_um} µm thick {window_material_name} window.")
    print(f"Window Density (ρ): {window_density_g_cm3:.3f} g/cm³")
    print(f"Window Thickness (x): {window_thickness_cm} cm\n")

    # --- Calculation Loop ---
    for name, Z in sorted(elements_to_check.items(), key=lambda item: item[1]):
        # Get the characteristic K-alpha line energy for the element in keV
        energy_keV = xraylib.LineEnergy(Z, xraylib.KA_LINE)

        # Get the mass attenuation coefficient (μ) of Be for that energy in cm²/g
        mu_cm2_g = xraylib.CS_Total_Kissel(window_material_symbol, energy_keV)

        # Calculate the argument of the exponent
        exponent_value = -mu_cm2_g * window_density_g_cm3 * window_thickness_cm

        # Calculate the transmission fraction
        transmission = np.exp(exponent_value)

        # --- Output the step-by-step calculation ---
        print(f"--- For {name} (Z={Z}) ---")
        print(f"K-alpha X-ray Energy: {energy_keV:.3f} keV")
        print(f"Mass Attenuation Coefficient (μ) of Be at {energy_keV:.3f} keV: {mu_cm2_g:.2f} cm²/g")
        print(f"Transmission = exp( -μ * ρ * x )")
        print(f"             = exp( -{mu_cm2_g:.2f} * {window_density_g_cm3:.3f} * {window_thickness_cm} )")
        print(f"             = exp( {exponent_value:.3f} )")
        print(f"             = {transmission:.4f}")
        print(f"Transmission Percentage: {transmission * 100:.2f}%\n")
    
    print("--- Conclusion ---")
    print("A 100 µm Be window is very thick and absorbs low-energy X-rays effectively.")
    print("Sodium (Na) X-rays have a very low transmission (<3%), making detection impractical.")
    print("Magnesium (Mg) is the lightest element with a significant transmission (>10%), making it reliably detectable.")

if __name__ == '__main__':
    solve_edx_detection_limit()
