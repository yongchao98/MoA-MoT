import math

def calculate_q_space():
    """
    Calculates the Q-space position for the second major diffraction peak
    of NaMgH3 perovskite.
    """
    # Step 1: Define the known parameters for NaMgH3 perovskite.
    # The crystal structure is simple cubic.
    # The lattice parameter 'a' at room temperature is approximately 3.99 Å.
    a = 3.99  # in Angstroms

    # Step 2: Identify the Miller indices (h, k, l) for the second major peak.
    # For a simple cubic perovskite, the first major peak is (110).
    # The second major peak observed in the diffraction pattern is (200).
    h, k, l = 2, 0, 0

    print("--- Calculation of Q-space Position ---")
    print(f"Material: NaMgH3 (perovskite)")
    print(f"Lattice parameter (a): {a} Å")
    print(f"Miller indices (h,k,l) for the second major peak: ({h}, {k}, {l})")
    print("-" * 35)

    # Step 3: Calculate the Q-space position using the formula.
    # Formula: Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)
    
    # Calculate the sum of the squares of the Miller indices
    hkl_sq_sum = h**2 + k**2 + l**2
    
    # Calculate the Q value
    q_value = (2 * math.pi / a) * math.sqrt(hkl_sq_sum)

    # Print the equation with the values substituted
    print("The Q-space position is calculated using the formula:")
    print("Q = (2 * \u03C0 / a) * \u221A(h\u00B2 + k\u00B2 + l\u00B2)")
    print("\nSubstituting the values:")
    print(f"Q = (2 * {math.pi:.4f} / {a}) * \u221A({h}\u00B2 + {k}\u00B2 + {l}\u00B2)")
    print(f"Q = ({2 * math.pi:.4f} / {a}) * \u221A({hkl_sq_sum})")
    
    print("\n--- Result ---")
    print(f"The Q-space position of the ({h},{k},{l}) peak is: {q_value:.4f} 1/Å")

# Execute the calculation
calculate_q_space()