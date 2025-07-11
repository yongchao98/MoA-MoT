import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3 perovskite.
    """
    # Step 1: Define crystal parameters for the pseudo-cubic NaMgH3 structure.
    # The lattice parameter 'a' in Angstroms.
    a = 4.13

    # Step 2: Identify Miller indices (h,k,l) for the second major diffraction peak.
    # For a perovskite, the first major peak is (110), the second is (200).
    h, k, l = 2, 0, 0
    
    print(f"Calculating Q-space position for the second major peak of NaMgH3, the ({h},{k},{l}) reflection.")
    print(f"Using a pseudo-cubic lattice parameter a = {a} Å.")
    print("-" * 30)

    # Step 3: Calculate the interplanar d-spacing.
    # Formula for a cubic system: d = a / sqrt(h^2 + k^2 + l^2)
    hkl_sum_sq = h**2 + k**2 + l**2
    d_spacing = a / math.sqrt(hkl_sum_sq)

    print(f"The d-spacing is calculated as: d = a / sqrt(h² + k² + l²)")
    print(f"d = {a} / sqrt({h}² + {k}² + {l}²) = {a} / sqrt({hkl_sum_sq}) = {d_spacing:.4f} Å")
    print("-" * 30)

    # Step 4: Calculate the Q-space position.
    # Formula: Q = 2 * pi / d
    Q = (2 * math.pi) / d_spacing

    print(f"The Q-space position is calculated as: Q = 2 * π / d")
    print(f"Q = (2 * {math.pi:.5f}) / {d_spacing:.4f} = {Q:.4f} 1/Å")

calculate_q_space_position()
