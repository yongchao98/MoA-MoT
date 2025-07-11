import math

def calculate_q_space_peak():
    """
    Calculates the Q-space position of the second major diffraction peak
    of NaMgH3 perovskite.
    """
    # Step 1: Define crystal parameters for NaMgH3
    # Lattice parameter 'a' for the idealized cubic structure in Angstroms.
    a = 3.89
    # Miller indices (h, k, l) for the second major peak in a simple cubic perovskite.
    h, k, l = 2, 0, 0

    print("This script calculates the Q-space position of the second major diffraction peak of NaMgH3.")
    print("------------------------------------------------------------------------------------------")
    print(f"1. Model: Ideal cubic perovskite structure.")
    print(f"2. Lattice parameter (a): {a} Å.")
    print(f"3. Miller indices (h,k,l) for the second major peak: ({h}, {k}, {l}).")
    print("------------------------------------------------------------------------------------------\n")

    # Step 2: Calculate the d-spacing for the (200) plane
    hkl_sq_sum = float(h**2 + k**2 + l**2)
    d_spacing = a / math.sqrt(hkl_sq_sum)

    # Step 3: Calculate the Q-space position
    q_value = (2 * math.pi) / d_spacing

    # Step 4: Print the calculation steps and final answer
    print("Calculation steps:")
    # Print d-spacing calculation
    print(f"d = a / sqrt(h² + k² + l²)")
    print(f"d = {a} / sqrt({h}² + {k}² + {l}²)")
    print(f"d = {a} / sqrt({hkl_sq_sum}) = {d_spacing:.4f} Å")
    print("")

    # Print Q-space calculation
    print(f"Q = 2 * π / d")
    print(f"Q = 2 * {math.pi:.4f} / {d_spacing:.4f}")
    print(f"Q = {q_value:.4f} Å⁻¹")
    print("\n------------------------------------------------------------------------------------------")
    print(f"The final Q-space position for the ({h},{k},{l}) peak is approximately {q_value:.2f} 1/Å.")


if __name__ == "__main__":
    calculate_q_space_peak()
    # The final numerical answer for the submission system.
    final_answer = 3.23 
    # print(f"<<<{final_answer}>>>") # This would be uncommented for submission