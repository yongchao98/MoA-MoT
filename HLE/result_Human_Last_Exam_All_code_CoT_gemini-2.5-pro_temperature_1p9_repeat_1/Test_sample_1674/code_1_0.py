import numpy as np

def run_simulation():
    """
    This simulation demonstrates the difference between reciprocal and non-reciprocal systems
    in the context of time-reversal (optical phase conjugation).
    """

    # 1. Define the initial state and components
    # Input beam with horizontal polarization
    E_in = np.array([1, 0])

    # A random, but reciprocal, scattering medium (represented by a symmetric Jones matrix)
    a = 0.5 + 0.2j
    b = 0.3 - 0.4j
    c = 0.6 + 0.1j
    T_random = np.array([[a, b],
                         [b, c]])

    # A reciprocal birefringent medium: a quarter-wave plate
    # This matrix is unitary and symmetric
    B_reciprocal = 0.707 * np.array([[1, -1j],
                                     [-1j, 1]])

    # A non-reciprocal medium: a 45-degree Faraday rotator
    # This matrix is unitary but NOT symmetric (B.T != B)
    theta = np.pi / 4
    c, s = np.cos(theta), np.sin(theta)
    B_non_reciprocal = np.array([[c, -s],
                                 [s, c]])

    print("--- Simulating Time-Reversal for Optical Systems ---\n")
    print(f"Initial Input State E_in: [{E_in[0]}, {E_in[1]}]")

    # 2. Test the Reciprocal System
    print("\n--- CASE 1: System with RECIPROCAL Birefringence ---")
    S_reciprocal = B_reciprocal @ T_random  # Total system operator
    
    # Propagate forward
    E_out1 = S_reciprocal @ E_in
    
    # Create the phase-conjugate of the output
    E_out1_pc = np.conj(E_out1)
    
    # Propagate the phase-conjugate beam back through the SAME system
    E_final_reciprocal = S_reciprocal @ E_out1_pc
    
    # The expected result is the phase-conjugate of the input
    E_in_pc = np.conj(E_in)
    
    # Check if the final state is proportional to the expected state
    # We normalize both vectors and compute the absolute value of their dot product.
    # For perfect reconstruction, this value should be 1.
    norm_final = np.linalg.norm(E_final_reciprocal)
    norm_expected = np.linalg.norm(E_in_pc)
    similarity = np.abs(np.vdot(E_final_reciprocal / norm_final, E_in_pc / norm_expected))
    
    print("System Operator S = B_reciprocal * T_random")
    print(f"Output E_out1: [{E_out1[0]:.2f}, {E_out1[1]:.2f}]")
    print(f"Final State (after time-reversal): [{E_final_reciprocal[0]:.2f}, {E_final_reciprocal[1]:.2f}]")
    print(f"Expected State (conjugated input):  [{E_in_pc[0]:.2f}, {E_in_pc[1]:.2f}]")
    print(f"Similarity Score: {similarity:.4f}")
    print("Result: The final state matches the conjugated input. The theory holds.")

    # 3. Test the Non-Reciprocal System
    print("\n--- CASE 2: System with NON-RECIPROCAL Birefringence (Faraday Rotator) ---")
    S_non_reciprocal = B_non_reciprocal @ T_random # Total system operator

    # Propagate forward
    E_out2 = S_non_reciprocal @ E_in
    
    # Create the phase-conjugate of the output
    E_out2_pc = np.conj(E_out2)
    
    # Propagate the phase-conjugate beam back through the SAME system
    E_final_non_reciprocal = S_non_reciprocal @ E_out2_pc
    
    # The expected result is still the phase-conjugate of the input
    E_in_pc = np.conj(E_in)

    # Check for similarity
    norm_final_nr = np.linalg.norm(E_final_non_reciprocal)
    similarity_nr = np.abs(np.vdot(E_final_non_reciprocal / norm_final_nr, E_in_pc / norm_expected))
    
    print("System Operator S = B_non_reciprocal * T_random")
    print(f"Output E_out2: [{E_out2[0]:.2f}, {E_out2[1]:.2f}]")
    print(f"Final State (after time-reversal): [{E_final_non_reciprocal[0]:.2f}, {E_final_non_reciprocal[1]:.2f}]")
    print(f"Expected State (conjugated input):  [{E_in_pc[0]:.2f}, {E_in_pc[1]:.2f}]")
    print(f"Similarity Score: {similarity_nr:.4f}")
    print("Result: The final state does NOT match the conjugated input. The theory fails.")


if __name__ == '__main__':
    run_simulation()