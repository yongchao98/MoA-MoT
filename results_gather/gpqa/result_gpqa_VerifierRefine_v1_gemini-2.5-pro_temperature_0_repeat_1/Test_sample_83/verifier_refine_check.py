import numpy as np
from scipy.signal import residue
from scipy.linalg import inv

def check_correctness_of_parallel_splitting_method():
    """
    This function verifies that the key factor for parallelizing the solution
    of du/dt = Au using a fractional approximation is the partial fraction
    decomposition of that approximation.

    It does this by:
    1. Setting up a sample problem (discretized 1D heat equation matrix A).
    2. Choosing a rational approximation R(z) for exp(z).
    3. Calculating the solution u(t+tau) in two ways:
        a) "Sequential": Directly as inv(Q(tau*A)) @ P(tau*A) @ u.
        b) "Parallel": By decomposing R(z) into partial fractions and solving
           the resulting independent linear systems.
    4. Comparing the results to show the methods are equivalent.
    5. Concluding that the partial fraction decomposition (Answer B) is the
       enabling mechanism.
    """

    # --- 1. Setup the Problem ---
    # Simulate the matrix 'A' from a 1D heat equation discretization.
    # The principle extends to higher dimensions.
    N = 10  # Number of spatial grid points
    dx = 1.0 / (N + 1)
    A = (np.diag(np.ones(N - 1), -1) - 2 * np.diag(np.ones(N), 0) + np.diag(np.ones(N - 1), 1)) / (dx**2)
    
    # Initial condition vector u(t=0)
    u0 = np.sin(np.pi * np.arange(1, N + 1) * dx)
    tau = 0.01  # Time step

    # --- 2. Define a Fractional Approximation R(z) = P(z)/Q(z) ---
    # We use the (1,2) Padé approximant to exp(z):
    # R(z) = (1 + z/3) / (1 - 2z/3 + z^2/6)
    # Numerator coefficients P(z) = p1*z + p0
    p_coeffs = [1/3, 1]
    # Denominator coefficients Q(z) = q2*z^2 + q1*z + q0
    q_coeffs = [1/6, -2/3, 1]

    # --- 3a. "Sequential" Calculation ---
    # This involves solving one large system: Q(tau*A) * u_new = P(tau*A) * u0
    identity = np.identity(N)
    P_matrix = p_coeffs[1] * identity + p_coeffs[0] * (tau * A)
    Q_matrix = q_coeffs[2] * identity + q_coeffs[1] * (tau * A) + q_coeffs[0] * (tau * A) @ (tau * A)
    
    # Solve for the new state u_sequential
    u_sequential = np.linalg.solve(Q_matrix, P_matrix @ u0)

    # --- 3b. "Parallel" Calculation using Partial Fractions (Answer B) ---
    # Decompose R(z) into: k + sum(r_j / (z - p_j))
    # Note: scipy.signal.residue needs coefficients of decreasing powers of z.
    try:
        residues, poles, k_term = residue(p_coeffs[::-1], q_coeffs[::-1])
    except Exception as e:
        return f"Could not perform partial fraction decomposition. Error: {e}"

    # The decomposition allows solving independent systems: (tau*A - p_j*I)x_j = u0
    # These solves can be done in parallel.
    solutions_x = []
    for pole in poles:
        # For each pole, form and solve a linear system.
        # Note: The poles and residues can be complex, so we use complex arithmetic.
        matrix_to_invert = tau * A - pole * identity
        x_j = np.linalg.solve(matrix_to_invert, u0)
        solutions_x.append(x_j)

    # Combine the results from the parallel solves
    u_parallel = k_term[0] * u0 if len(k_term) > 0 else 0
    for i in range(len(residues)):
        u_parallel += residues[i] * solutions_x[i]

    # --- 4. Compare Results ---
    # The results should be numerically identical (within a small tolerance).
    if not np.allclose(u_sequential, u_parallel):
        return (f"Incorrect. The parallel decomposition method did not yield the same result as the sequential method. "
                f"This would mean that the principle described in answer B is flawed. "
                f"Max difference: {np.max(np.abs(u_sequential - u_parallel))}")

    # --- 5. Conclusion ---
    # The successful demonstration confirms that the partial fraction decomposition is the
    # mathematical mechanism that converts the problem into independent, parallelizable tasks.
    # - A) Stability is a necessary property of the chosen R(z), but not the mechanism for parallelism.
    # - C) Nonlocal boundary conditions would complicate matrix A, but not change the decomposition principle.
    # - D) The nature of the roots (complex in this case) is an implementation detail; the decomposition itself is the key.
    
    return "Correct"

# Execute the check
result = check_correctness_of_parallel_splitting_method()

# Final check based on the code's conclusion
if result == "Correct":
    print("Correct")
else:
    print(f"Incorrect: {result}")
