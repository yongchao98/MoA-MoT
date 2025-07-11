import numpy as np

def random_walk_analysis(N=8):
    """
    Analyzes a 1D random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    if N < 2:
        print("N must be an integer greater than or equal to 2.")
        return

    print(f"--- Analysis of a Random Walk on a Circle with N = {N} sites ---")

    # Part 1: One-step transformation equation
    print("\n1. One-Step Transformation Equation:")
    print("Let P_t(i) be the probability at site i. The distribution at time t+1 is:")
    print(f"P_t+1(i) = 0.5 * P_t((i-1) mod N) + 0.5 * P_t((i+1) mod N)")

    # Part 2: Transition Probability Matrix A
    print("\n2. Transition Probability Matrix (A):")
    # A[i, j] is the probability of transition from state j to i
    A = np.zeros((N, N))
    for j in range(N):
        # Move right: j -> (j+1) % N
        A[(j + 1) % N, j] = 0.5
        # Move left: j -> (j-1) % N
        A[(j - 1) % N, j] = 0.5
    print("The A matrix is:")
    print(A)

    # Part 3: Verification of Eigenvalues
    print("\n3. Eigenvalue Verification:")
    print("Analytical formula: λ_n = cos(2 * π * n / N) for n = 0, ..., N-1")
    # Calculate analytical eigenvalues and sort them
    analytical_eigs = np.sort([np.cos(2 * np.pi * n / N) for n in range(N)])[::-1]
    
    # Calculate numerical eigenvalues and sort them
    numerical_eigs, _ = np.linalg.eig(A)
    # Take real part to discard negligible imaginary parts from numerical errors
    numerical_eigs = np.sort(numerical_eigs.real)[::-1]

    print("\nComparison of numerically computed vs. analytical eigenvalues (sorted):")
    print("-" * 45)
    print(f"{'Numerical':>20} | {'Analytical':>20}")
    print("-" * 45)
    for num_eig, an_eig in zip(numerical_eigs, analytical_eigs):
        print(f"{num_eig:20.10f} | {an_eig:20.10f}")
    print("-" * 45)
    if np.allclose(numerical_eigs, analytical_eigs):
        print("Verification successful: The eigenvalues match.")
    else:
        print("Verification failed: Eigenvalues do not match.")

    # Part 4: Rate of Relaxation
    print("\n4. Rate of Relaxation:")
    print("The rate is determined by the second-largest eigenvalue (λ_1).")
    second_largest_eig_val = analytical_eigs[1]
    pi_val = np.pi
    
    print("\nThe equation for the second-largest eigenvalue λ_1 is:")
    print(f"λ_1 = cos(2 * π / N)")
    print(f"    = cos(2 * {pi_val:.6f} / {N})")
    print(f"    = cos({(2 * pi_val / N):.6f})")
    print(f"    = {second_largest_eig_val:.15f}")
    print("\n--------------------------------------------------------------")
    
    # Return the final value for the specified format
    return second_largest_eig_val

# Run the analysis for N=8
final_rate = random_walk_analysis(N=8)

# The final answer in the requested format
print(f"\nFinal Answer: The rate of relaxation determined by the second-largest eigenvalue is {final_rate}")
<<<0.7071067811865476>>>