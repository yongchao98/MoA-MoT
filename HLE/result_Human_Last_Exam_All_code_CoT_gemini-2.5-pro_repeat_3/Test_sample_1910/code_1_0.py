import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak
    of NaMgH3 perovskite.
    """
    # Step 1: Define lattice parameters from a reliable source (Materials Project).
    # Units are in Angstroms (Å).
    a = 5.859
    b = 8.277
    c = 5.920

    # Step 2: Define the Miller indices (h, k, l) for the second major peak.
    # Based on simulated XRD patterns, this is the (202) peak.
    h, k, l = 2, 0, 2

    print("Calculating the Q-space position for the (hkl) = (202) peak of orthorhombic NaMgH3.")
    print(f"Lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å.")
    print("The formula is: Q = sqrt((h/a)^2 + (k/b)^2 + (l/c)^2)\n")

    # Step 3: Calculate each term in the equation.
    h_a_sq = (h / a)**2
    k_b_sq = (k / b)**2
    l_c_sq = (l / c)**2
    
    print("Calculation steps:")
    print(f"(h/a)^2 = ({h}/{a})^2 = {h_a_sq:.6f}")
    print(f"(k/b)^2 = ({k}/{b})^2 = {k_b_sq:.6f}")
    print(f"(l/c)^2 = ({l}/{c})^2 = {l_c_sq:.6f}")

    # Step 4: Sum the terms to get Q^2.
    q_squared = h_a_sq + k_b_sq + l_c_sq
    print(f"\nQ^2 = {h_a_sq:.6f} + {k_b_sq:.6f} + {l_c_sq:.6f} = {q_squared:.6f} Å^-2")

    # Step 5: Take the square root to find Q.
    q_value = math.sqrt(q_squared)
    print(f"Q = sqrt({q_squared:.6f}) = {q_value:.4f} Å^-1\n")
    
    print(f"The Q-space position of the second major diffraction peak (202) is located at {q_value:.4f} 1/Å.")
    
    # Returning the final value for the answer tag
    return q_value

if __name__ == '__main__':
    final_q_value = calculate_q_space_position()
    # The final answer is formatted below as requested.
    # print(f"<<<{final_q_value:.4f}>>>")