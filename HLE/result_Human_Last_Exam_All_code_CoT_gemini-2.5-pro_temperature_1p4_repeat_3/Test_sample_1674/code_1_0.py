import numpy as np

def run_simulation():
    """
    This simulation demonstrates how adding a birefringent medium can cause
    the user's described theory to fail.
    """
    # --- Step 1: Define the optical components and input ---

    # Set a random seed for reproducibility
    np.random.seed(0)

    # 1a. An arbitrary input beam with horizontal polarization
    # [Ex, Ey] where Ex is the x-component and Ey is the y-component.
    E_in = np.array([1, 0], dtype=complex)

    # 1b. A random, invertible 2x2 matrix representing the random medium (T)
    # This matrix can scramble the polarization.
    T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
    T_inv = np.linalg.inv(T)

    # 1c. A Jones matrix for a birefringent medium (B)
    # This represents a quarter-wave plate with its fast axis at 45 degrees.
    # B^2 will be the matrix for a half-wave plate, so B^2 != B.
    B = 0.5 * np.array([[1-1j, 1+1j], [1+1j, 1-1j]], dtype=complex)


    print("--- Simulation Setup ---")
    print(f"Input Polarization Vector (E_in):\n{E_in}\n")
    print(f"Random Medium Matrix (T):\n{np.round(T, 2)}\n")
    print(f"Birefringent Plate Matrix (B):\n{np.round(B, 2)}\n")


    # --- Step 2: Test the theory WITHOUT the birefringent plate ---
    # The system is just the random medium T.

    print("--- Case 1: System without Birefringence (Theory Holds) ---")

    # The output from the system
    E_out1 = T @ E_in
    print(f"Original Output (E_out1 = T @ E_in):\n{np.round(E_out1, 2)}\n")

    # Create the "necessary input" using the user's method
    E_in2 = T_inv @ E_out1
    # Note: E_in2 should be identical to E_in, proving the method finds the original input
    # print(f"Calculated Input (E_in2 = T_inv @ E_out1):\n{np.round(E_in2, 2)}\n")


    # Test the theory: does feeding E_in2 back into the system yield E_out1?
    Test_Output = T @ E_in2
    print("Equation Check: T @ E_in2 == E_out1")
    print(f"{np.round(Test_Output, 2)} == {np.round(E_out1, 2)}")
    print(f"Are they equal? {np.allclose(Test_Output, E_out1)}\n\n")


    # --- Step 3: Test the theory WITH the birefringent plate ---
    # The new system operator is M = B @ T

    print("--- Case 2: System with Birefringence (Theory Fails) ---")

    # The new system operator
    M = B @ T

    # The output from the new system
    E_out1_new = M @ E_in
    print(f"New Original Output (E_out1_new = B @ T @ E_in):\n{np.round(E_out1_new, 2)}\n")

    # Create the "necessary input" using the same method as before
    # (inverting only the T part of the system)
    E_in2_new = T_inv @ E_out1_new

    # Test the theory: does feeding E_in2_new back into the NEW system (M) yield E_out1_new?
    Test_Output_new = M @ E_in2_new

    print("Final Equation Check: M @ E_in2_new == E_out1_new")
    print("This is equivalent to checking if: (B^2 @ T @ E_in) == (B @ T @ E_in)\n")
    print("Calculated final result (LHS):")
    print(np.round(Test_Output_new, 2))
    print("\nOriginal desired result (RHS):")
    print(np.round(E_out1_new, 2))

    print(f"\nAre they equal? {np.allclose(Test_Output_new, E_out1_new)}")


if __name__ == '__main__':
    run_simulation()