import sympy as sp

def analyze_random_walk_on_circle():
    """
    Analyzes the 1D random walk on a circle to find its transition matrix,
    eigenvalues, and the rate of relaxation.
    """
    # Use unicode for prettier printing in compatible terminals
    sp.init_printing(use_unicode=True)

    # --- 1. One-step Transformation ---
    print("--- 1. One-step Transformation ---")
    print("Consider a random walk on a circle with N sites, indexed j = 0, 1, ..., N-1.")
    print("Let p_j(t) be the probability of being at site j at time t.")
    print("At each step, the walker moves to a neighboring site with probability 1/2.")
    print("The one-step transformation of the probability distribution is given by the equation:")

    p_j_t1 = sp.Symbol('p_j(t+1)')
    p_j_minus_1_t = sp.Symbol('p_{j-1}(t)')
    p_j_plus_1_t = sp.Symbol('p_{j+1}(t)')
    transformation_eq = sp.Eq(p_j_t1, sp.S(1)/2 * p_j_minus_1_t + sp.S(1)/2 * p_j_plus_1_t)
    print("\nFinal Equation 1: Transformation Rule")
    sp.pprint(transformation_eq)
    print("\nThis means the probability at site j at the next step is the average of the probabilities at its neighbors.")

    # --- 2. Transition Matrix A ---
    print("\n\n--- 2. The Transition Matrix A ---")
    print("This transformation can be written in matrix form as p(t+1) = A * p(t).")
    print("The matrix element A_ij = P(move to site i | at site j).")
    print("Therefore, A_ij = 1/2 if i = (j-1)%N or i = (j+1)%N, and 0 otherwise.")
    print("For example, with N=4, the transition matrix A is:")
    A_example = sp.Matrix([
        [0, sp.S(1)/2, 0, sp.S(1)/2],
        [sp.S(1)/2, 0, sp.S(1)/2, 0],
        [0, sp.S(1)/2, 0, sp.S(1)/2],
        [sp.S(1)/2, 0, sp.S(1)/2, 0]
    ])
    sp.pprint(A_example)

    # --- 3. Eigenvector and Eigenvalue Calculation ---
    print("\n\n--- 3. Eigenvalue and Eigenvector Derivation ---")
    N = sp.Symbol('N', integer=True, positive=True)
    j = sp.Symbol('j', integer=True)
    n = sp.Symbol('n', integer=True)

    print("We solve the eigenvalue equation A * v_n = lambda_n * v_n.")
    print("The proposed eigenvectors v_n have components v_nj = exp(i * k_n * j),")
    print(f"where the wave number k_n = 2*pi*n/N for n = 0, 1, ..., N-1.\n")

    k_n = 2 * sp.pi * n / N
    v_nj = sp.exp(sp.I * k_n * j)

    print("Let's apply the transformation A to the j-th component of v_n:")
    print("(A * v_n)_j = (1/2) * v_{n, j-1} + (1/2) * v_{n, j+1}")

    v_n_j_minus_1 = sp.exp(sp.I * k_n * (j - 1))
    v_n_j_plus_1 = sp.exp(sp.I * k_n * (j + 1))
    result_j = sp.S(1)/2 * (v_n_j_minus_1 + v_n_j_plus_1)

    print(f"            = (1/2) * (exp(i*k_n*(j-1)) + exp(i*k_n*(j+1)))")
    print(f"            = exp(i*k_n*j) * (1/2) * (exp(-i*k_n) + exp(i*k_n))")

    # Using Euler's formula, cos(x) = (e^(ix) + e^(-ix))/2
    eigenvalue_expr = sp.cos(k_n)
    print(f"            = exp(i*k_n*j) * cos(k_n)")
    print(f"            = {sp.pretty(eigenvalue_expr)} * {sp.pretty(v_nj)}\n")
    
    print("This confirms that v_nj are the components of the eigenvectors.")
    print("The corresponding eigenvalues lambda_n are:")
    lambda_n_sym = sp.Symbol('lambda_n')
    lambda_eq = sp.Eq(lambda_n_sym, eigenvalue_expr)

    print("\nFinal Equation 2: Eigenvalues")
    sp.pprint(lambda_eq)

    # --- 4. Rate of Relaxation ---
    print("\n\n--- 4. Rate of Relaxation ---")
    print("The system relaxes towards the stationary distribution, which corresponds to the largest eigenvalue, lambda_0 = 1.")
    lambda_0 = eigenvalue_expr.subs(n, 0)
    print(f"For n = 0, lambda_0 = cos(0) = {lambda_0}")

    print("\nThe rate of relaxation is determined by the slowest decaying mode, which corresponds to the second-largest eigenvalue.")
    print("This is the eigenvalue with the largest magnitude less than 1.")
    print("The eigenvalues are cos(2*pi*n/N). These values are largest when n is small.")
    print("The second-largest distinct eigenvalue occurs for n=1 (and n=N-1):")
    
    lambda_1 = eigenvalue_expr.subs(n, 1)
    lambda_2nd_largest = sp.Symbol('lambda_relax')
    relax_eq = sp.Eq(lambda_2nd_largest, lambda_1)

    print("\nFinal Equation 3: Eigenvalue Determining Relaxation Rate")
    sp.pprint(relax_eq)


if __name__ == '__main__':
    analyze_random_walk_on_circle()
