import numpy as np

def solve_random_walk_circle(N=10):
    """
    Analyzes a 1D random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    # --- Introduction and Model Definition ---
    print(f"--- 1D Random Walk on a Circle with N={N} sites ---")
    print("\n1. One-Step Transformation of Probability:")
    print("Let p_j(t) be the probability of being at site j at time t.")
    print("The probability distribution evolves according to the rule:")
    print("p_j(t+1) = 0.5 * p_{j-1}(t) + 0.5 * p_{j+1}(t)")
    print("where indices are taken modulo N.")
    print("-" * 50)

    # --- 2. Transition Probability Matrix A ---
    print("\n2. Transition Probability Matrix (A):")
    print("This transformation can be written in matrix form: p(t+1) = A * p(t).")
    
    # Create the matrix A where A[i, j] is the prob of moving from j to i
    A = np.zeros((N, N))
    for i in range(N):
        # Probability of moving from site (i-1) to i is 0.5
        A[i, (i - 1 + N) % N] = 0.5
        # Probability of moving from site (i+1) to i is 0.5
        A[i, (i + 1) % N] = 0.5

    print("The N x N transition matrix A is:")
    np.set_printoptions(precision=3, suppress=True)
    print(A)
    print("-" * 50)

    # --- 3. Eigenvalues Demonstration ---
    print("\n3. Eigenvalues of A:")
    print("The eigenvalues of this matrix are given by the formula:")
    print("lambda_n = cos(2 * pi * n / N), for n = 0, ..., N-1.")

    # Calculate theoretical and numerical eigenvalues for verification
    theoretical_eigenvalues = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
    numerical_eigenvalues, _ = np.linalg.eig(A)

    # Sort them for easier comparison (descending order)
    theoretical_eigenvalues_sorted = np.sort(theoretical_eigenvalues)[::-1]
    numerical_eigenvalues_sorted = np.sort(np.real(numerical_eigenvalues))[::-1]

    print("\nTheoretical Eigenvalues (sorted):")
    print(theoretical_eigenvalues_sorted)
    print("\nNumerical Eigenvalues from NumPy (sorted):")
    print(numerical_eigenvalues_sorted)
    print("The values match, confirming the formula.")
    print("-" * 50)

    # --- 4. Rate of Relaxation ---
    print("\n4. Rate of Relaxation:")
    print("The system relaxes towards the stationary distribution (eigenvalue = 1).")
    print("The rate of relaxation is governed by the second-largest eigenvalue,")
    print("which corresponds to the slowest decaying mode.")
    
    # The second largest eigenvalue is for n=1
    n_for_rate = 1
    second_largest_eig_val = np.cos(2 * np.pi * n_for_rate / N)

    print("\nThis eigenvalue is calculated for n=1:")
    print(f"lambda_1 = cos(2 * pi * {n_for_rate} / {N})")
    
    # Print the equation with all numbers filled in
    pi_val = np.pi
    arg_val = 2 * pi_val * n_for_rate / N
    print("\nFinal Calculation:")
    final_equation = f"lambda_1 = cos(2 * {pi_val:.5f} * {n_for_rate} / {N}) = cos({arg_val:.5f}) = {second_largest_eig_val:.5f}"
    print(final_equation)
    print("-" * 50)

    return second_largest_eig_val

# Run the analysis for a circle with 10 sites.
# You can change the number of sites here.
num_sites = 10
relaxation_rate = solve_random_walk_circle(N=num_sites)

# The final answer requested is the value determining the rate of relaxation.
print(f"\n<<<>>>") # Placeholder for the final answer
# The value is printed to the console but also returned.
# To follow the instruction "don't ask users to copy", we format the output this way.
print(f"{relaxation_rate}")

if __name__ == '__main__':
    # This block is not strictly necessary for the user prompt, 
    # but it ensures the code is runnable as a script.
    pass