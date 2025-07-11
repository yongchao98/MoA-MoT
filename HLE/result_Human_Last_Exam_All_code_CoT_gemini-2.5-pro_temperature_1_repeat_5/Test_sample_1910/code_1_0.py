import math

def calculate_q_space_peak():
    """
    Calculates the Q-space position of the second major diffraction peak of NaMgH3.
    """
    # Step 1: Define lattice parameters for orthorhombic NaMgH3 (Pnma) in Angstroms.
    a = 5.46
    b = 7.74
    c = 5.39

    # Step 2: Identify the Miller indices (h, k, l) for the second major peak.
    # The first major peak is a cluster ((011), (101), etc.).
    # The second major peak is commonly the (022) reflection.
    h, k, l = 0, 2, 2

    print("Calculating the Q-space position for the (hkl) = ({}, {}, {}) peak of NaMgH3.".format(h, k, l))
    print("Lattice parameters (orthorhombic): a = {} Å, b = {} Å, c = {} Å\n".format(a, b, c))

    # Step 3: Calculate the interplanar spacing (d_hkl).
    # The formula for an orthorhombic system is: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2
    h_a_sq = (h / a)**2
    k_b_sq = (k / b)**2
    l_c_sq = (l / c)**2
    d_inv_sq = h_a_sq + k_b_sq + l_c_sq
    
    # Avoid division by zero if d_inv_sq is zero
    if d_inv_sq == 0:
        print("Error: Miller indices are (0,0,0), which is not a diffraction peak.")
        return

    d_hkl = 1 / math.sqrt(d_inv_sq)

    # Step 4: Calculate the Q value.
    # The formula is: Q = 2 * pi / d_hkl
    Q = 2 * math.pi / d_hkl

    # Step 5: Print the detailed calculation.
    print("The formula for interplanar spacing 'd' is: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2")
    print("1/d^2 = ({}/{})^2 + ({}/{})^2 + ({}/{})^2".format(h, a, k, b, l, c))
    print("1/d^2 = {:.4f} + {:.4f} + {:.4f} = {:.4f} Å⁻²".format(h_a_sq, k_b_sq, l_c_sq, d_inv_sq))
    print("d = 1 / sqrt({:.4f}) = {:.4f} Å\n".format(d_inv_sq, d_hkl))

    print("The formula for the scattering vector 'Q' is: Q = 2 * π / d")
    print("Q = 2 * π / {:.4f} Å".format(d_hkl))
    print("Q = {:.4f} Å⁻¹".format(Q))
    
    # Return the final Q value for the answer block
    return Q

if __name__ == '__main__':
    q_value = calculate_q_space_peak()
    # The final answer is formatted as requested by the user.
    # print(f"\n<<<{q_value:.3f}>>>")
