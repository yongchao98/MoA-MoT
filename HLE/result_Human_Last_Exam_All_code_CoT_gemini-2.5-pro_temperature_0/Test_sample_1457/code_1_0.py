import math

def solve_edx_detection():
    """
    Calculates the lightest element detectable by an EDX system with a thick Be window.
    """
    # --- Physical Constants ---
    # 1. Candidate elements with their atomic number (Z) and K-alpha X-ray energy
    elements = {
        "Na": {"Z": 11, "Ka_energy_keV": 1.041},
        "Mg": {"Z": 12, "Ka_energy_keV": 1.254},
        "Si": {"Z": 14, "Ka_energy_keV": 1.740},
        "Ca": {"Z": 20, "Ka_energy_keV": 3.692},
    }

    # 2. Beryllium (Be) window properties
    be_density_g_cm3 = 1.85
    window_thickness_cm = 100e-4  # 100 µm = 0.01 cm

    # 3. Mass attenuation coefficients (μ/ρ) for Be at the elements' Kα energies (from NIST XCOM database)
    # Units: cm²/g
    be_attenuation_coeffs = {
        1.041: 1233.0,  # For Na Kα
        1.254: 743.5,   # For Mg Kα
        1.740: 333.1,   # For Si Kα
        3.692: 47.6,    # For Ca Kα
    }

    # --- Calculation ---
    print("Calculating X-ray transmission through a 100 µm Beryllium window...")
    print("-" * 70)

    detectable_elements = []
    # Sort elements by atomic number to find the lightest one first
    sorted_elements = sorted(elements.items(), key=lambda item: item[1]['Z'])

    for name, properties in sorted_elements:
        energy = properties["Ka_energy_keV"]
        z_number = properties["Z"]
        mu_rho = be_attenuation_coeffs[energy]

        # Beer-Lambert Law: T = exp(- (μ/ρ) * ρ * x)
        exponent_value = -mu_rho * be_density_g_cm3 * window_thickness_cm
        transmission = math.exp(exponent_value)

        print(f"Element: {name} (Z={z_number})")
        print(f"Kα Energy: {energy} keV")
        print(f"Calculation: T = exp( -(μ/ρ) * ρ * x )")
        print(f"             T = exp( -({mu_rho} cm²/g) * ({be_density_g_cm3} g/cm³) * ({window_thickness_cm} cm) )")
        print(f"             T = exp( {exponent_value:.4f} )")
        print(f"Transmission: {transmission:.6f} (or {transmission * 100:.4f}%)")

        # Check against the 1% detection threshold
        if transmission > 0.01:
            detectable_elements.append(name)
            print("Result: DETECTABLE (Transmission > 1%)\n")
        else:
            print("Result: NOT DETECTABLE (Transmission <= 1%)\n")

    # --- Conclusion ---
    print("-" * 70)
    if detectable_elements:
        lightest_element = detectable_elements[0]
        print(f"The lightest element with transmission > 1% is: {lightest_element}")
    else:
        print("None of the candidate elements are detectable with a 1% threshold.")

solve_edx_detection()
<<<D>>>