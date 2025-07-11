import numpy as np

def random_walk_on_circle(N):
    """
    Analyzes the random walk on a circle with N sites.

    This function performs the following steps:
    1. Defines the one-step probability transformation.
    2. Constructs the transition probability matrix A.
    3. Shows that the eigenvectors are of the form exp(i*j*k_n) and finds the eigenvalues.
    4. Identifies the second-largest eigenvalue which determines the relaxation rate.
    """
    print(f"--- Analysis for a 1D Random Walk on a Circle with N = {N} sites ---")
    print("-" * 60)

    # 1. One-Step Transformation
    print("Step 1: The One-Step Transformation\n")
    print("Let p(t) be the probability distribution vector at time t.")
    print("The distribution at the next step, p(t+1), is given by a matrix multiplication:")
    print("p(t+1) = A * p(t)\n")
    print("where A is the N x N transition probability matrix.\n")
    print("-" * 60)

    # 2. Transition Probability Matrix A
    print("Step 2: Constructing the Transition Matrix A\n")
    print("A[j, i] is the probability of moving from site i to site j.")
    print("From any site i, the walker can move to i-1 or i+1 (modulo N) with probability 1/2 each.")
    
    A = np.zeros((N, N))
    for i in range(N):
        A[(i - 1) % N, i] = 0.5
        A[(i + 1) % N, i] = 0.5
    
    print("\nThe transition matrix A for N = {} is:\n".format(N))
    print(A)
    print("\n" + "-" * 60)

    # 3. Eigenvectors and Eigenvalues
    print("Step 3: Verifying Eigenvectors and Deriving Eigenvalues\n")
    print("We want to show that vectors v_n with components v_n(j) = exp(i*j*k_n)")
    print(f"where i is the imaginary unit, k_n = 2*pi*n/N, and n = 0, 1, ..., N-1, are eigenvectors of A.\n")
    print("We test the eigenvector equation: A * v_n = lambda_n * v_n")
    print("Let's check the j-th component of the left side, (A * v_n)_j:")
    print("\n(A * v_n)_j = Sum over l [ A[j, l] * v_n(l) ]")
    print(f"           = A[j, j-1]*v_n(j-1) + A[j, j+1]*v_n(j+1)   (indices are mod {N})")
    print("           = 0.5 * v_n(j-1) + 0.5 * v_n(j+1)")
    print("\nSubstitute v_n(j) = exp(i*j*k_n):")
    print("           = 0.5 * exp(i*(j-1)*k_n) + 0.5 * exp(i*(j+1)*k_n)")
    print("Factor out exp(i*j*k_n):")
    print("           = 0.5 * exp(i*j*k_n) * [ exp(-i*k_n) + exp(i*k_n) ]")
    print("Using Euler's formula, exp(ix) + exp(-ix) = 2*cos(x):")
    print("           = 0.5 * exp(i*j*k_n) * [ 2*cos(k_n) ]")
    print("           = cos(k_n) * exp(i*j*k_n)")
    print("           = cos(k_n) * v_n(j)")
    print("\nSo, A * v_n = cos(k_n) * v_n. This is true for all components j.")
    print("This confirms that v_n are eigenvectors, and the eigenvalues are:")
    print("lambda_n = cos(k_n) = cos(2*pi*n / N) for n = 0, 1, ..., N-1.\n")
    print("-" * 60)
    
    # 4. Relaxation Rate
    print("Step 4: Finding the Rate of Relaxation\n")
    print("The rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    print("The largest eigenvalue corresponds to n=0:")
    print("lambda_0 = cos(2*pi*0 / N) = cos(0) = 1.0")
    print("This corresponds to the stationary distribution (a uniform distribution).\n")
    print("The second-largest eigenvalue corresponds to n=1 (and n=N-1), because cos(x) is largest when x is close to 0.")
    print("Therefore, the determining factor for the relaxation rate is lambda_1.\n")
    
    # Final Equation and Calculation
    print("Final Equation for the second-largest eigenvalue:")
    second_largest_eigenvalue_formula = "cos(2 * pi / N)"
    n_val = 1
    second_largest_eigenvalue_val = np.cos(2 * np.pi * n_val / N)

    print(f"lambda_1 = {second_largest_eigenvalue_formula}\n")

    print(f"For our specific case with N = {N}:")
    print(f"lambda_1 = cos(2 * pi / {N}) = {second_largest_eigenvalue_val:.8f}")

if __name__ == '__main__':
    # You can change the number of sites N here
    N = 5
    random_walk_on_circle(N)
