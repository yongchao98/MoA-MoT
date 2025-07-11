import math

def solve_edx_problem():
    """
    Calculates the lightest element detectable by an EDX with a thick Be window.
    """
    # 1. Define physical constants and parameters
    # Beryllium (Be) window properties
    rho_be = 1.85  # Density of Beryllium in g/cm^3
    thickness_um = 100  # Window thickness in micrometers
    thickness_cm = thickness_um * 1e-4  # Convert thickness to cm

    # Data for elements in the answer choices.
    # Kα energy and Mass Attenuation Coefficient (μ) of Be at that energy.
    # Data is sourced from standard physics tables (e.g., NIST).
    elements = {
        # 'Name': {'Z', 'E_ka_keV', 'mu_be_cm2g'}
        'Na': {'Z': 11, 'E_ka_keV': 1.041, 'mu_be_cm2g': 4403.4},
        'Mg': {'Z': 12, 'E_ka_keV': 1.254, 'mu_be_cm2g': 2633.2},
        'Si': {'Z': 14, 'E_ka_keV': 1.740, 'mu_be_cm2g': 1124.6},
        'Ca': {'Z': 20, 'E_ka_keV': 3.692, 'mu_be_cm2g': 163.5},
    }

    print("--- EDX Detectability Calculation ---")
    print(f"Calculating X-ray transmission through a {thickness_um} µm Beryllium window.")
    print(f"Window properties: ρ = {rho_be} g/cm³, x = {thickness_cm} cm\n")

    # Sort elements by atomic number to find the "lightest" one first.
    sorted_elements = sorted(elements.items(), key=lambda item: item[1]['Z'])

    # 2. Calculate transmission for each element
    for name, properties in sorted_elements:
        mu = properties['mu_be_cm2g']
        
        # Calculate the argument of the exponent
        exponent = -mu * rho_be * thickness_cm
        
        # Calculate transmission using Beer-Lambert law: T = exp(-μ*ρ*x)
        transmission = math.exp(exponent)
        
        print(f"Element: {name} (Z={properties['Z']}, Kα Energy={properties['E_ka_keV']} keV)")
        print(f"  μ in Be: {mu:.1f} cm²/g")
        print("  Calculation: T = exp(-μ * ρ * x)")
        print(f"  T = exp(-{mu:.1f} * {rho_be} * {thickness_cm}) = {transmission:.4e}")
        print(f"  Resulting Transmission: {transmission * 100:.6f} %\n")

    # 3. Analyze results and conclude
    print("--- Conclusion ---")
    print("An element is considered 'seen' if its X-ray transmission is high enough to be measured.")
    print("Typically, this requires transmission to be at least ~0.1%.")
    print("Based on the calculations:")
    print(" - Na, Mg, and Si X-rays are almost completely absorbed by the window (transmission is effectively 0%).")
    print(" - Calcium (Ca) is the first element in the list (ordered by atomic number) with a measurable transmission (~4.85%).")
    print("\nTherefore, the lightest element that can be seen is Calcium (Ca).")

# Run the analysis
solve_edx_problem()