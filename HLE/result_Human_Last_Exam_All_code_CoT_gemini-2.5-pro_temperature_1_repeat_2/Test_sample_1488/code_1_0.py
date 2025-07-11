import numpy as np

def analyze_xrd_pattern():
    """
    Analyzes XRD data to determine phase, structure, and lattice constants.
    This function uses pre-identified peaks from the provided data to calculate
    the lattice parameters for a hypothesized tetragonal structure.
    """
    # --- Input Parameters ---
    # Wavelength of the X-ray source in Angstroms (standard for Cu K-alpha)
    wavelength = 1.5406

    # Prominent peaks identified from the XRD data (in degrees)
    # These correspond to the (101) and (200) planes of Rutile TiO2
    peak1_2theta = 35.98
    peak2_2theta = 39.18
    
    # Miller indices (h,k,l) for the identified peaks
    (h1, k1, l1) = (1, 0, 1)
    (h2, k2, l2) = (2, 0, 0)

    # --- Analysis ---
    print("--- XRD Analysis ---")
    print("Based on the prominent peaks, the material is identified as the Rutile phase of Titanium Dioxide (TiO2).")
    print("\nIdentification Details:")
    print(f"Assumed X-ray Wavelength (Cu Kα): {wavelength} Å")
    print("Hypothesized Crystal System: Tetragonal")
    print(f"Peak 1 at 2θ = {peak1_2theta}° is indexed as ({h1}, {k1}, {l1}).")
    print(f"Peak 2 at 2θ = {peak2_2theta}° is indexed as ({h2}, {k2}, {l2}).")

    # --- Calculations ---
    # Relationship for a tetragonal system: 1/d² = (h² + k²)/a² + l²/c²
    
    # Convert angles from degrees to radians
    theta1_rad = np.deg2rad(peak1_2theta / 2.0)
    theta2_rad = np.deg2rad(peak2_2theta / 2.0)

    # Calculate d-spacing using Bragg's Law: d = λ / (2 * sin(θ))
    d1 = wavelength / (2 * np.sin(theta1_rad))
    d2 = wavelength / (2 * np.sin(theta2_rad))

    # Calculate 1/d^2 values
    one_over_d1_sq = 1 / (d1**2)
    one_over_d2_sq = 1 / (d2**2)

    # Solve for lattice constants 'a' and 'c'
    # From Peak 2 (200): 1/d₂² = (2² + 0²)/a² + 0²/c²  =>  1/d₂² = 4/a²
    a_sq = (h2**2 + k2**2) / one_over_d2_sq
    a = np.sqrt(a_sq)

    # From Peak 1 (101): 1/d₁² = (1² + 0²)/a² + 1²/c²  =>  1/d₁² = 1/a² + 1/c²
    # Rearranging for c: 1/c² = 1/d₁² - 1/a²
    one_over_c_sq = one_over_d1_sq - (1 / a_sq)
    c_sq = 1 / one_over_c_sq
    c = np.sqrt(c_sq)

    print("\n--- Calculation Steps ---")
    print("Step 1: Calculate lattice constant 'a' from the (2 0 0) peak.")
    print(f"Equation: 1 / d² = ( {h2}² + {k2}² ) / a²")
    print(f"For 2θ = {peak2_2theta:.2f}°, d = {d2:.4f} Å, so 1/d² = {one_over_d2_sq:.4f} Å⁻².")
    print(f"a² = ( {h2}² + {k2}² ) / {one_over_d2_sq:.4f} = {a_sq:.4f} Å²")
    print(f"a = √{a_sq:.4f} = {a:.3f} Å")
    
    print("\nStep 2: Calculate lattice constant 'c' using 'a' and the (1 0 1) peak.")
    print(f"Equation: 1 / c² = (1 / d²) - ( ( {h1}² + {k1}² ) / a² )")
    print(f"For 2θ = {peak1_2theta:.2f}°, d = {d1:.4f} Å, so 1/d² = {one_over_d1_sq:.4f} Å⁻².")
    print(f"1 / c² = {one_over_d1_sq:.4f} - ( ( {h1}² + {k1}² ) / {a_sq:.4f} ) = {one_over_c_sq:.4f} Å⁻²")
    print(f"c² = 1 / {one_over_c_sq:.4f} = {c_sq:.4f} Å²")
    print(f"c = √{c_sq:.4f} = {c:.3f} Å")

    print("\n\n--- Final Results ---")
    print(f"Chemical Phase: Titanium Dioxide (TiO₂), Rutile structure")
    print(f"Unit Cell Structure: Tetragonal")
    print(f"Lattice Constants: a = {a:.3f} Å, c = {c:.3f} Å")

if __name__ == '__main__':
    analyze_xrd_pattern()
<<<Chemical Phase: Titanium Dioxide (TiO₂), Rutile structure
Unit Cell Structure: Tetragonal
Lattice Constants: a = 4.593 Å, c = 2.971 Å>>>