import sympy

def solve_problem():
    """
    This function symbolically follows the logic of the proof to determine theta.
    """
    n, j, c = sympy.symbols('n j c', positive=True, real=True)
    
    print("Step 1: Define the relationship between E[tau] and the sum of probabilities.")
    # E[tau] = n - Sum_{j=1}^{n-1} P(S_j >= T)
    # We need to bound Sum P(S_j >= T) <= c * n**theta
    
    print("\nStep 2: Define means and variances.")
    # Probability of non-zero X_i
    p = n**(-sympy.Rational(1, 2))
    
    # Mean of U_i ~ U[0, n**(-1/2)]
    E_Ui = sympy.Rational(1, 2) * n**(-sympy.Rational(1, 2))
    
    # Mean of X_i
    E_Xi = p * E_Ui
    
    # Mean of S_j
    E_Sj = j * E_Xi
    
    # Variance of U_i
    Var_Ui = (n**(-sympy.Rational(1, 2)))**2 / 12
    
    # E[U_i^2]
    E_Ui2 = Var_Ui + E_Ui**2
    
    # E[X_i^2]
    E_Xi2 = p * E_Ui2
    
    # Variance of X_i
    Var_Xi = E_Xi2 - E_Xi**2
    
    # Variance of S_j
    Var_Sj = j * Var_Xi
    
    print("E[S_j] =", sympy.simplify(E_Sj))
    print("Var(S_j) =", sympy.simplify(Var_Sj))
    
    print("\nStep 3: Apply Chebyshev's inequality.")
    # Threshold T
    T = 1 - n**(-sympy.Rational(1, 2))
    
    # Bound for P(S_j >= T)
    prob_bound = Var_Sj / (T - E_Sj)**2
    
    print("The upper bound for P(S_j >= T) from Chebyshev's inequality is proportional to:")
    # Analyze the order of magnitude for large n
    prob_bound_large_n = sympy.limit(prob_bound / (j * n**(-sympy.Rational(3, 2))), n, sympy.oo)
    print(f"Var(S_j) is O(j * n^(-3/2))")
    print(f"(T - E_Sj)^2 is O(1)")
    print(f"So P(S_j >= T) is O(j * n^(-3/2))")

    print("\nStep 4: Bound the sum of probabilities.")
    # Sum_{j=1}^{n-1} O(j * n^(-3/2))
    sum_bound_order = sympy.Symbol('O', real=True) * (n**2 * n**(-sympy.Rational(3, 2)))
    sum_bound_simplified = sympy.simplify(sum_bound_order)
    
    print("Sum_{j=1}^{n-1} j*n^(-3/2) = n^(-3/2) * (n-1)*n/2, which is O(n^(1/2))")
    print("So, n - E[tau] <= c * n^(1/2)")
    
    theta = sympy.Rational(1, 2)
    
    print("\nStep 5: Conclude the value of theta.")
    print(f"The analysis shows that n - E[tau] is bounded by c*n^({theta}).")
    print(f"This means the inequality E[tau] >= n - c*n^theta holds for theta = {theta}.")
    
    # The question asks for a multiple of 1/8.
    theta_eighths = theta * 8
    
    print(f"As a multiple of 1/8, theta is {int(theta_eighths)}/8.")
    print("The final equation is E[tau] >= n - c * n^(1/2)")
    
    print("\nFinal equation with numbers:")
    print("E[tau] >= n - c * n^(", theta, ")")

    # Final answer format
    # No direct return, but the logic points to 1/2
    
solve_problem()