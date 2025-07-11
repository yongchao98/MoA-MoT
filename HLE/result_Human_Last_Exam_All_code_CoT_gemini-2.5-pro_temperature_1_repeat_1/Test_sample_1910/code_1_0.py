import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3.
    """
    # Step 1: Define lattice parameters for orthorhombic NaMgH3 at room temperature (in Angstroms)
    a = 5.518
    b = 7.794
    c = 5.509

    # Step 2: Define Miller indices (h, k, l) for the second major peak.
    # The first major peak is (111). The second is a cluster of (121), (200), and (002).
    # We will use (1, 2, 1) as the representative indices for this peak.
    h, k, l = 1, 2, 1
    
    print(f"Crystal: NaMgH3 (orthorhombic)")
    print(f"Lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å")
    print(f"Miller indices for the second major peak: (h, k, l) = ({h}, {k}, {l})\n")

    # Step 3: Calculate the d-spacing
    # Using the formula: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2
    h_a_sq = (h / a)**2
    k_b_sq = (k / b)**2
    l_c_sq = (l / c)**2
    d_inv_sq = h_a_sq + k_b_sq + l_c_sq
    d_spacing = math.sqrt(1 / d_inv_sq)
    
    print("--- Calculation Steps ---")
    print("1. Calculate d-spacing:")
    print(f"1/d² = (h²/a²) + (k²/b²) + (l²/c²)")
    print(f"1/d² = ({h}²/{a}²) + ({k}²/{b}²) + ({l}²/{c}²)")
    print(f"1/d² = {h_a_sq:.5f} + {k_b_sq:.5f} + {l_c_sq:.5f} = {d_inv_sq:.5f} Å⁻²")
    print(f"d = 1 / sqrt({d_inv_sq:.5f}) = {d_spacing:.4f} Å\n")

    # Step 4: Calculate the Q-space vector magnitude
    # Using the formula: Q = 2 * pi / d
    q_value = 2 * math.pi / d_spacing

    print("2. Calculate Q-space position:")
    print(f"Q = 2 * π / d")
    print(f"Q = 2 * {math.pi:.5f} / {d_spacing:.4f}")
    print(f"Q = {q_value:.4f} Å⁻¹\n")

    print("--- Final Answer ---")
    print(f"The second major diffraction peak is located at a Q-space position of {q_value:.4f} 1/Å.")
    
    # Return the final numeric answer for the platform
    global final_answer
    final_answer = round(q_value, 4)


# Run the calculation
calculate_q_space_position()
<<<2.2792>>>