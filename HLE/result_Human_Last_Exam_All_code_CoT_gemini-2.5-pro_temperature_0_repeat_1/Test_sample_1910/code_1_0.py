import math

def calculate_diffraction_peak_q_space():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3.

    The crystal structure is orthorhombic (Pnma) with lattice parameters:
    a = 5.53 Å
    b = 7.80 Å
    c = 5.51 Å

    The second major peak (ordered by increasing Q-value) corresponds to the
    (101) and (020) Miller indices. We will perform the calculation for (101).
    """
    # Lattice parameters in Angstroms
    a = 5.53
    b = 7.80
    c = 5.51

    # Miller indices for the second peak
    h, k, l = 1, 0, 1

    print("Calculating the Q-space position for the second major diffraction peak of NaMgH3.")
    print(f"The peak corresponds to Miller indices (h,k,l) = ({h},{k},{l}).")
    print(f"Lattice parameters (orthorhombic): a = {a} Å, b = {b} Å, c = {c} Å\n")

    # Calculate terms for the 1/d^2 equation
    h_a_sq = (h / a)**2
    k_b_sq = (k / b)**2
    l_c_sq = (l / c)**2
    d_inv_sq = h_a_sq + k_b_sq + l_c_sq

    # Calculate Q
    q_value = 2 * math.pi * math.sqrt(d_inv_sq)

    # Print the step-by-step calculation
    print("The formula for the inverse squared d-spacing is:")
    print("1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2")
    print(f"1/d^2 = ({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2")
    print(f"1/d^2 = {h_a_sq:.5f} + {k_b_sq:.5f} + {l_c_sq:.5f}")
    print(f"1/d^2 = {d_inv_sq:.5f} Å⁻²\n")

    print("The formula for the Q-space vector is:")
    print("Q = 2 * pi * sqrt(1/d^2)")
    print(f"Q = 2 * {math.pi:.5f} * sqrt({d_inv_sq:.5f})")
    print(f"Q = {q_value:.3f} Å⁻¹")

calculate_diffraction_peak_q_space()
<<<1.610>>>