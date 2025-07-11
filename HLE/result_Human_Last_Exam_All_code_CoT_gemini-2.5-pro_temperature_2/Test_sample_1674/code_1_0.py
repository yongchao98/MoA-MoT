import numpy as np

def main():
    """
    Demonstrates that adding a birefringent element to an optical system
    invalidates a pre-calculated inverse of the original system.
    """
    # Set a random seed for reproducibility
    np.random.seed(0)
    np.set_printoptions(precision=3, suppress=True)

    # --- 1. Define the Initial Components ---

    # Define the input beam as horizontally polarized light.
    # E = [Ex, Ey]
    E_in = np.array([[1 + 0j], [0 + 0j]])
    print("--- System Setup ---")
    print(f"Initial Input Beam (Horizontal Polarization):\n{E_in}\n")

    # Create a random 2x2 Jones matrix for the "random medium".
    # This matrix scrambles the polarization state of the input beam.
    T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
    print(f"Random Medium's Jones Matrix (T):\n{T}\n")

    # Calculate the inverse of this matrix. This is our "magic" correction plate.
    try:
        T_inv = np.linalg.inv(T)
        print(f"Inverse of the Random Medium (T_inv):\n{T_inv}\n")
    except np.linalg.LinAlgError:
        print("The random matrix T is singular and cannot be inverted. Please run again.")
        return

    # --- 2. The Original Theory: Inverting the Random Medium ---

    print("\n--- Part 1: Verifying the Original Theory (No Birefringence) ---")
    # Propagate the beam through the random medium
    E_out_original = T @ E_in
    print(f"Output after passing through T:\n{E_out_original}\n")

    # Apply the inverse matrix to recover the original beam
    E_recovered_original = T_inv @ E_out_original
    print(f"Result of applying T_inv to the output:")
    print(f"Recovered Beam:\n{E_recovered_original}\n")

    # Check if the recovered beam is the same as the input
    if np.allclose(E_in, E_recovered_original):
        print("Success! The original beam was recovered. The theory holds for the original system.\n")
    else:
        print("Error in calculation. The original beam was not recovered.\n")


    # --- 3. The New System: Adding a Birefringent Plate ---

    print("\n--- Part 2: Testing the Theory after Adding a Birefringent Plate ---")
    # Define a Jones matrix for a birefringent medium.
    # This is a quarter-wave plate, which induces a pi/2 phase shift between
    # the x and y components.
    # Let's assume the light passes through this plate *before* the random medium.
    J_B = np.array([[1, 0], [0, 1j]])
    print(f"Birefringent Plate's Jones Matrix (B):\n{J_B}\n")

    # Propagate the initial beam through the NEW system (Plate B, then Medium T)
    E_after_B = J_B @ E_in
    print(f"Beam state after plate B (but before T):\n{E_after_B}\n")

    E_out_new_system = T @ E_after_B
    print(f"Output after the full new system (B then T):\n{E_out_new_system}\n")


    # --- 4. Attempt to Correct the New System with the OLD Inverse ---

    print("Attempting to correct the new output using the OLD inverse (T_inv)...")
    # Apply the original inverse T_inv to the output of the new system
    E_recovered_badly = T_inv @ E_out_new_system
    print(f"Final Recovered Beam:\n{E_recovered_badly}\n")

    # Check if the recovered beam is the same as the original input
    if np.allclose(E_in, E_recovered_badly):
        print("This is unexpected. The original beam was recovered.")
    else:
        print("Failure! The recovered beam does NOT match the original input beam.")
        print("The theory does not hold for the new system.")

    # We can see that the badly recovered beam is equal to the state *after* the
    # birefringent plate. The T_inv only undid the effect of T.
    assert np.allclose(E_recovered_badly, E_after_B)
    print("\nNote: The final state is equal to the state of the beam just after the birefringent plate. The old inverse T_inv successfully undid the effect of T, but the effect of B was left uncorrected.")

if __name__ == "__main__":
    main()
