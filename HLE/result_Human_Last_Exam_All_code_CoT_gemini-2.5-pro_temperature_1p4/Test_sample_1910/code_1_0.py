import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3.
    """
    # Step 1: Define the pseudo-cubic lattice parameter for NaMgH3 in Angstroms.
    # This is derived from the experimental orthorhombic cell parameters.
    a_p = 3.895

    # Step 2: Define the Miller indices (h,k,l) for the second major peak.
    # For a perovskite, this is typically the (200) reflection.
    h, k, l = 2, 0, 0

    # Step 3: Calculate the d-spacing for the (h,k,l) peak.
    d_spacing = a_p / math.sqrt(h**2 + k**2 + l**2)

    # Step 4: Calculate the position in Q-space.
    Q = (2 * math.pi) / d_spacing

    # Print the explanation and the calculation step-by-step.
    print(f"Based on the pseudo-cubic model for NaMgH3 perovskite:")
    print(f"Lattice parameter 'a' = {a_p} Å")
    print(f"Miller indices (h,k,l) for the second major peak = ({h},{k},{l})")
    print("\nCalculation of d-spacing:")
    print(f"d = a / sqrt(h^2 + k^2 + l^2)")
    print(f"d = {a_p} / sqrt({h}^2 + {k}^2 + {l}^2) = {d_spacing:.4f} Å")
    
    print("\nCalculation of Q-space position:")
    print(f"Q = 2 * pi / d")
    print(f"Q = 2 * {math.pi:.4f} / {d_spacing:.4f} = {Q:.3f} 1/Å")
    
    # Return the final numeric value for the answer tag.
    return Q

# Run the calculation and store the result.
q_value = calculate_q_space_position()

# The final answer will be extracted from this line by the system.
# print(f"<<<{q_value:.3f}>>>")