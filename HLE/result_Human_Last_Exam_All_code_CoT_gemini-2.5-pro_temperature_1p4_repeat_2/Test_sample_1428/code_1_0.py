import numpy as np

def analyze_random_walk_on_circle(N):
    """
    Constructs the transition matrix for a 1D random walk on a circle,
    computes its eigenvalues, and finds the rate of relaxation.

    Args:
        N (int): The number of sites on the circle.
    """
    if not isinstance(N, int) or N <= 2:
        print("Error: N must be an integer greater than 2.")
        return

    print(f"--- Analysis for Random Walk on a Circle with N = {N} sites ---")

    # 1. Construct the Transition Matrix A
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5
    
    print("\n[Step 1] The transition matrix A is:")
    print(A)

    # 2. Compute Eigenvalues
    # Numerical computation
    numerical_eigvals = np.linalg.eigvals(A)
    # Sort in descending order for easier comparison
    numerical_eigvals_sorted = np.sort(numerical_eigvals)[::-1]

    # Analytical computation
    analytical_eigvals = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
    analytical_eigvals_sorted = np.sort(analytical_eigvals)[::-1]
    
    print("\n[Step 2] Verifying the eigenvalues (sorted in descending order):")
    print(f"Numerical:   {[f'{v:.5f}' for v in numerical_eigvals_sorted]}")
    print(f"Analytical:  {[f'{v:.5f}' for v in analytical_eigvals_sorted]}")
    assert np.allclose(numerical_eigvals_sorted, analytical_eigvals_sorted), "Mismatch found!"
    print("Verification successful: Numerical and analytical eigenvalues match.")

    # 3. Find the Rate of Relaxation
    # The largest eigenvalue corresponds to n=0
    lambda_0 = analytical_eigvals_sorted[0]
    
    # The second-largest eigenvalue corresponds to n=1
    second_largest_eigval = analytical_eigvals_sorted[1]

    print("\n[Step 3] Finding the rate of relaxation:")
    print("The rate is determined by the second-largest eigenvalue (λ_1).")
    print(f"Largest eigenvalue (λ_0 for n=0): cos(2*π*0/{N}) = {lambda_0:.5f}")
    
    print("\nThe final equation for the second-largest eigenvalue is:")
    # Using '1' for n=1 explicitly in the equation string
    final_equation_str = f"λ_1 = cos(2 * π * 1 / {N})"
    final_value = second_largest_eigval
    print(f"{final_equation_str} = {final_value:.6f}")


if __name__ == '__main__':
    # Define the number of sites for the random walk.
    # You can change this value to see the results for different circle sizes.
    number_of_sites = 10
    analyze_random_walk_on_circle(number_of_sites)
