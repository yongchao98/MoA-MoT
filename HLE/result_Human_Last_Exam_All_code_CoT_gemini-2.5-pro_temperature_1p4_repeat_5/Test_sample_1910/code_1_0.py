import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak
    of NaMgH3 perovskite.
    """
    # Step 1: Define the constants for the calculation.
    # Lattice parameter 'a' for NaMgH3 in Angstroms.
    a = 3.90
    # Miller indices (h, k, l) for the second major diffraction peak ({200}).
    h, k, l = 2, 0, 0

    print("--- Calculation of Q-space Position ---")
    print(f"Material: NaMgH3 (cubic perovskite)")
    print(f"Lattice parameter (a): {a} Å")
    print(f"Miller indices (h,k,l) for the second major peak: ({h},{k},{l})\n")

    # Step 2: Calculate the Q value using the formula.
    # Formula: Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)
    hkl_squared_sum = h**2 + k**2 + l**2
    q_value = (2 * math.pi / a) * math.sqrt(hkl_squared_sum)

    # Step 3: Print the detailed calculation and the final result.
    print("Formula: Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)")
    print("Substituting the values:")
    # Printing each number in the final equation as requested
    print(f"Q = (2 * {math.pi} / {a}) * sqrt({h}^2 + {k}^2 + {l}^2)")
    print(f"Q = {2 * math.pi / a} * sqrt({hkl_squared_sum})")
    print(f"Q = {2 * math.pi / a} * {math.sqrt(hkl_squared_sum)}")
    print(f"Final calculated Q-space position: {q_value:.4f} Å⁻¹")

if __name__ == "__main__":
    calculate_q_space_position()
