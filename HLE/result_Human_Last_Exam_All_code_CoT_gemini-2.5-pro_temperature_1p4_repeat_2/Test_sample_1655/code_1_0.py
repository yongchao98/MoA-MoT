import numpy as np

def demonstrate_vector_beam_limitation():
    """
    Numerically demonstrates that an N-DoF input cannot generate an arbitrary 2N-DoF output.
    """
    # 1. Define the number of spatial modes (degrees of freedom, DoF) for the input.
    # This represents the number of pixels or modes we can independently control.
    N = 10

    print(f"Input Degrees of Freedom (DoF): {N}")
    print(f"Required Output Degrees of Freedom (DoF) for an arbitrary vector beam: {2 * N}\n")

    # 2. Model the optical system with a random Transmission Matrix (T).
    # This matrix T maps an input scalar field (N complex values) to an output vector
    # field (2N complex values, for Ex and Ey at each of N points).
    # The shape of T is (2N, N).
    np.random.seed(0)  # Use a seed for reproducible results
    T = np.random.randn(2 * N, N) + 1j * np.random.randn(2 * N, N)
    print("Created a random Transmission Matrix T of shape (2N, N) to model the optical system.")

    # 3. Create a desired arbitrary "target" vector beam.
    # This is a random vector in the 2N-dimensional output space that we want to create.
    E_target = np.random.randn(2 * N) + 1j * np.random.randn(2 * N)
    print("Created a random 'target' vector beam E_target.\n")

    # 4. Attempt to find the input S_in required to generate E_target.
    # We need to solve the equation: T @ S_in = E_target.
    # This is a system of 2N equations for N unknowns, which is an overdetermined system.
    # We can find the least-squares solution, which finds the S_in that minimizes the
    # error || T @ S_in - E_target ||^2.
    print("Solving for the input S_in that best produces E_target...")
    # np.linalg.lstsq finds the least-squares solution.
    S_in_solution, residuals, rank, s = np.linalg.lstsq(T, E_target, rcond=None)

    # 5. Calculate the actual beam produced and the error.
    E_produced = T @ S_in_solution
    error_norm_squared = np.sum(np.abs(E_target - E_produced)**2)

    print("\n--- Results ---")
    print(f"Norm-squared of target beam ||E_target||^2 = {np.sum(np.abs(E_target)**2):.4f}")

    # Outputting the error calculation in an equation-like format
    print(f"\nFinal Equation for the error:")
    print(f"||T * S_in - E_target||^2 = {error_norm_squared:.4f}")

    # 6. Conclusion
    # A non-zero error means the target beam could not be perfectly generated.
    if not np.isclose(error_norm_squared, 0):
        print("\nThe error is significantly non-zero.")
        print("This demonstrates that the arbitrary target vector beam is not reachable with the given system.")
    else:
        # This case is statistically almost impossible for random matrices and targets.
        print("\nThe error is zero. The arbitrary target was reachable (highly unlikely).")

if __name__ == '__main__':
    demonstrate_vector_beam_limitation()