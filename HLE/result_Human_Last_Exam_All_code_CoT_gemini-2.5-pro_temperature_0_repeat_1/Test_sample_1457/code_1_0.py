import math
try:
    import xraylib
except ImportError:
    print("Please install the xraylib library using: pip install xraylib")
    exit()

def solve_edx_detection_limit():
    """
    Calculates the transmission of characteristic X-rays from light elements
    through a Beryllium detector window to find the lightest detectable element.
    """
    # --- Problem Parameters ---
    # Elements to check from the answer choices (Name, Symbol, Atomic Number Z)
    # W is the sample, not a light element impurity we are looking for.
    elements_to_check = [
        ("Sodium", "Na", 11),
        ("Magnesium", "Mg", 12),
        ("Silicon", "Si", 14),
        ("Calcium", "Ca", 20)
    ]

    # Beryllium window parameters
    be_atomic_number = 4
    window_thickness_um = 100.0
    window_thickness_cm = window_thickness_um / 10000.0  # Convert µm to cm

    # Get physical constants for Beryllium using xraylib
    be_density_g_cm3 = xraylib.ElementDensity(be_atomic_number)

    print(f"Analyzing X-ray transmission through a {window_thickness_um} µm Be window (Density: {be_density_g_cm3:.2f} g/cm³)\n")
    print("-" * 70)
    print(f"{'Element':<12} | {'Z':<3} | {'Kα Energy (keV)':<16} | {'Transmission (%)':<18}")
    print("-" * 70)

    results = []

    for name, symbol, z in elements_to_check:
        # Get the K-alpha line energy for the element in keV
        # xraylib.LineEnergy returns energy in keV
        ka_energy_kev = xraylib.LineEnergy(z, xraylib.KA_LINE)

        # Get the mass attenuation coefficient of Be for this energy in cm²/g
        # xraylib.CS_MassCoef takes Z of absorber and energy in keV
        mass_atten_coeff_cm2_g = xraylib.CS_MassCoef(be_atomic_number, ka_energy_kev)

        # Calculate transmission using the Beer-Lambert law: I/I₀ = exp(-μ * ρ * x)
        exponent = -mass_atten_coeff_cm2_g * be_density_g_cm3 * window_thickness_cm
        transmission = math.exp(exponent)
        transmission_percent = transmission * 100

        results.append({
            "symbol": symbol,
            "transmission": transmission,
            "mass_atten_coeff": mass_atten_coeff_cm2_g
        })

        print(f"{name:<12} | {z:<3} | {ka_energy_kev:<16.3f} | {transmission_percent:<18.6f}")

    print("-" * 70)
    print("\n--- Conclusion ---")
    print("Sodium (Na) and Magnesium (Mg) X-rays are almost completely absorbed by the window.")
    print("Silicon (Si) is the lightest element in the list with a non-negligible (though small) transmission.")
    print("Therefore, Si is the lightest element that can be seen.\n")

    # Find the result for Silicon to print the final equation
    si_result = next(item for item in results if item["symbol"] == "Si")
    si_transmission = si_result["transmission"]
    si_mu = si_result["mass_atten_coeff"]

    print("Final calculation for the lightest detectable element (Si):")
    # Final equation with all numbers printed
    print(f"Transmission = exp(-μ * ρ * x)")
    print(f"Transmission = exp(-{si_mu:.2f} cm²/g * {be_density_g_cm3:.2f} g/cm³ * {window_thickness_cm:.4f} cm)")
    print(f"Transmission = {si_transmission*100:.4f}%")


solve_edx_detection_limit()
<<<E>>>