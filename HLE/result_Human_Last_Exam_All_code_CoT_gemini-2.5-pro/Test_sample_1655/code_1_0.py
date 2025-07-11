import numpy as np

def demonstrate_dof_limitation():
    """
    Demonstrates the mathematical limitation of generating an arbitrary vector
    beam from a scalar input beam.
    """
    # Let's model the beam's cross-section with N discrete points (pixels).
    N = 10

    # 1. Input Space (Scalar Beams)
    # The input is a scalar beam, with controlled amplitude and phase at each point.
    # This corresponds to N complex degrees of freedom.
    # The dimension of the input space is N.
    input_dim = N

    # 2. Output Space (Arbitrary Vector Beams)
    # An arbitrary vector beam requires two independent complex numbers (e.g.,
    # for horizontal and vertical polarization) at each point.
    # This corresponds to 2*N complex degrees of freedom.
    # The dimension of the full output space is 2*N.
    output_dim = 2 * N

    print("--- Analysis of Degrees of Freedom (DoF) ---")
    print(f"Number of spatial points (N): {input_dim}")
    print(f"Input DoF (Scalar Beam): {input_dim}")
    print(f"Required Output DoF (Arbitrary Vector Beam): {output_dim}")
    print("-" * 45)

    # 3. The System as a Linear Transformation
    # The optical system, no matter how complex, acts as a linear operator.
    # It maps the N-dimensional input space to the 2N-dimensional output space.
    # We can model this transformation with a random matrix of shape (2*N, N).
    # This matrix represents the entire system: propagations, random medium, etc.
    system_matrix = np.random.randn(output_dim, input_dim) + 1j * np.random.randn(output_dim, input_dim)

    # 4. Dimensionality of Achievable Outputs
    # The set of all possible output beams is the column space of the system matrix.
    # The dimension of this space is the rank of the matrix. For a (2N, N) matrix
    # with random entries, the rank will be N.
    achievable_dim = np.linalg.matrix_rank(system_matrix)

    print("\n--- System Capability Analysis ---")
    print(f"The optical system maps the {input_dim}-DoF input to a {output_dim}-DoF space.")
    print(f"The dimension of all achievable outputs (rank of system matrix) is: {achievable_dim}")

    if achievable_dim < output_dim:
        print("\nConclusion: Since the dimension of achievable outputs is smaller than the")
        print("dimension of the target space, it is IMPOSSIBLE to generate an ARBITRARY")
        print("vector beam. We can only generate beams within a limited subspace.")
    else:
        # This case is mathematically impossible here.
        print("\nConclusion: The system can generate any arbitrary vector beam.")

    # 5. Numerical Demonstration
    # Let's create an arbitrary "target" vector beam we want to generate.
    # This is just a random vector in the 2N-dimensional space.
    target_vector_beam = np.random.randn(output_dim) + 1j * np.random.randn(output_dim)

    # Now, let's try to find an input scalar beam that produces this target.
    # We are trying to solve the equation: system_matrix * input = target
    # This is an overdetermined system (more equations than unknowns).
    # We use a least-squares solver to find the input that gives the closest
    # possible output to our target.
    _, residuals, _, _ = np.linalg.lstsq(system_matrix, target_vector_beam, rcond=None)

    print("\n--- Numerical Demonstration ---")
    print("Attempting to create an arbitrary target vector beam...")
    # The residual measures the squared error between the best possible output and the target.
    # A non-zero residual means the target cannot be perfectly created.
    # Since residuals is an array, we access the first (and only) element.
    error = residuals[0]
    print(f"Least-squares error (residual) for a random target: {error:.4f}")
    if error > 1e-9:
        print("The non-zero error proves that the arbitrary target cannot be generated.")
    else:
        print("The target was generated (this is statistically almost impossible).")

if __name__ == '__main__':
    demonstrate_dof_limitation()
    print("\nBased on this analysis, the final answer is No.")