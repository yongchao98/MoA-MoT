import math

def calculate_q_space_peak():
    """
    Calculates the Q-space position for the second major diffraction peak of orthorhombic NaMgH3.
    """
    # Step 1: Define lattice parameters for orthorhombic NaMgH3 (in Angstroms)
    a = 5.922
    b = 8.324
    c = 5.869

    # Step 2: Define Miller indices (h, k, l) for the second major peak
    h, k, l = 2, 1, 1

    print(f"Calculating the Q-space position for the ({h},{k},{l}) peak of NaMgH3.")
    print(f"Using lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å\n")

    # Step 3: Calculate the interplanar spacing (d)
    # 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2
    h_a_sq = (h / a)**2
    k_b_sq = (k / b)**2
    l_c_sq = (l / c)**2
    d_sq_inv = h_a_sq + k_b_sq + l_c_sq
    
    # Avoid division by zero, although not an issue here
    if d_sq_inv == 0:
        print("Error: Miller indices are (0,0,0) or lattice parameters are invalid.")
        return

    d = math.sqrt(1 / d_sq_inv)

    print("Step 1: Calculate d-spacing")
    print(f"1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2")
    print(f"1/d^2 = ({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2")
    print(f"1/d^2 = {h_a_sq:.5f} + {k_b_sq:.5f} + {l_c_sq:.5f} = {d_sq_inv:.5f} Å⁻²")
    print(f"d = sqrt(1 / {d_sq_inv:.5f}) = {d:.4f} Å\n")
    
    # Step 4: Calculate the Q-space position
    # Q = 2 * pi / d
    q_value = 2 * math.pi / d
    
    print("Step 2: Calculate Q-space position")
    print(f"Q = 2 * π / d")
    print(f"Q = 2 * {math.pi:.4f} / {d:.4f} Å = {q_value:.3f} Å⁻¹")


calculate_q_space_peak()
print("<<<2.494>>>")