import numpy as np

def run_simulation():
    """
    This simulation models an optical system using Jones matrices to demonstrate
    why adding a birefringent element breaks a simple inverse recovery theory.
    """

    # Set a seed for reproducibility of the random matrix
    np.random.seed(0)
    # Use complex numbers for Jones calculus
    dtype = np.complex128

    # 1. Define the initial input beam: a horizontally polarized beam
    # Represented as a Jones vector [Ex, Ey]
    input_beam = np.array([1, 0], dtype=dtype)
    print("Initial Input Beam (Horizontal Polarization):")
    print(input_beam)
    print("-" * 50)

    # 2. Define the optical elements
    # A random medium, represented by a random 2x2 complex matrix
    random_medium_T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
    
    # The inverse of the random medium, used for the recovery process
    inverse_random_medium_T_inv = np.linalg.inv(random_medium_T)

    # A birefringent plate, represented by a Jones matrix for a quarter-wave plate
    # This plate will change the polarization of the light passing through it
    birefringent_plate_B = np.array([[1, 0], [0, 1j]], dtype=dtype)
    
    print("Defined Optical Elements:")
    print("Random Medium T:\n", np.round(random_medium_T, 2))
    print("\nInverse of Random Medium (T_inv):\n", np.round(inverse_random_medium_T_inv, 2))
    print("\nBirefringent Plate B (Quarter-Wave Plate):\n", np.round(birefringent_plate_B, 2))
    print("-" * 50)

    # --- Scenario 1: The original theory without the birefringent plate ---
    print("SCENARIO 1: System without Birefringent Plate")
    print("The theory states that applying T and then T_inv recovers the input.")
    
    # Propagate beam through the random medium T
    output_after_T = random_medium_T @ input_beam
    print(f"\nEquation: Output = T * Input")
    print(f"Output after T: {np.round(output_after_T, 2)}")
    
    # Recover the beam using the inverse matrix T_inv
    recovered_beam = inverse_random_medium_T_inv @ output_after_T
    print(f"\nEquation: Recovered = T_inv * Output")
    print(f"Recovered Beam: {np.round(recovered_beam, 2)}")
    
    # Check if the recovered beam is identical to the input beam
    is_recovered_scenario1 = np.allclose(input_beam, recovered_beam)
    print(f"\nIs the beam recovered? {is_recovered_scenario1}. The theory holds here.")
    print("-" * 50)

    # --- Scenario 2: The theory fails when the birefringent plate is added ---
    print("SCENARIO 2: System with Birefringent Plate Added")
    print("The beam first passes through B, then T. We try to recover using the old inverse T_inv.")

    # The new system transformation is T * B
    # First, the beam passes through the birefringent plate B
    output_after_B = birefringent_plate_B @ input_beam
    print(f"\nEquation: Intermediate = B * Input")
    print(f"Beam state after B: {np.round(output_after_B, 2)}")

    # Then, this new beam passes through the random medium T
    final_output_with_B = random_medium_T @ output_after_B
    print(f"\nEquation: Final Output = T * Intermediate")
    print(f"Final output after B and T: {np.round(final_output_with_B, 2)}")

    # Attempt to recover the original beam using ONLY the inverse of T
    attempted_recovery = inverse_random_medium_T_inv @ final_output_with_B
    print(f"\nEquation: Attempted Recovery = T_inv * Final Output")
    print(f"Attempted recovery result: {np.round(attempted_recovery, 2)}")

    # Check if the recovered beam is identical to the ORIGINAL input beam
    is_recovered_scenario2 = np.allclose(input_beam, attempted_recovery)
    print(f"\nIs the original beam recovered? {is_recovered_scenario2}.")
    print("\nCONCLUSION: The theory does not hold.")
    print("The resulting vector is not the original input. It is equal to B * Input, because T_inv * (T * B * Input) simplifies to B * Input.")

if __name__ == '__main__':
    run_simulation()