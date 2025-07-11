import sympy

def solve_random_walk_problem():
    """
    Solves the given random walk problem by interpreting the escape probability,
    finding its asymptotic behavior, and then calculating the limit.
    """

    # Define the symbols for the symbolic calculation
    n = sympy.Symbol('n', real=True, positive=True)
    c = sympy.Symbol('c', real=True, positive=True) # Proportionality constant

    # Step 1 & 2: Explain the interpretation and find the scaling of p_n
    print("Step 1: Interpretation of the escape probability p_n")
    print("The problem considers a simple random walk on the 3D integer lattice Z^3.")
    print("A 3D random walk is transient, meaning it eventually wanders off to infinity.")
    print("If 'p_n' were the probability of just leaving the cube C_n, p_n would be 1, making the limit 0.")
    print("A more meaningful interpretation is that p_n is the probability that the walk leaves C_n and *never returns*.")
    print("\nStep 2: Asymptotic behavior of p_n")
    print("The probability of escaping to infinity from a large object of size L in 3D is proportional to 1/L.")
    print("The cube C_n has a characteristic size L that is proportional to n.")
    print("Therefore, p_n asymptotically behaves as c/n for some constant c.")
    
    p_n_asymptotic = c / n
    print(f"\nAsymptotic form: p_n ≈ {p_n_asymptotic}")

    # Step 3: Set up and evaluate the limit
    print("\nStep 3: Calculating the limit of ln(1/p_n) / ln(n)")
    
    # The expression to be evaluated
    expression = sympy.log(1 / p_n_asymptotic) / sympy.log(n)
    
    print("The expression is: lim_{n->inf} [ ln(1 / (c/n)) / ln(n) ]")
    
    # Show simplification steps which makes the final equation clear
    print("\nSimplifying the expression gives us the final equation steps:")
    # Using specific numbers for the equation as requested.
    # We follow the symbolic simplification of the limit.
    print("lim_{n->inf} (ln(n/c) / ln(n))")
    print("= lim_{n->inf} ( (ln(n) - ln(c)) / ln(n) )")
    print("= lim_{n->inf} ( 1 - ln(c)/ln(n) )")
    
    # As n -> inf, ln(n) -> inf, so ln(c)/ln(n) -> 0
    # The final numbers in the equation are 1 and 0.
    final_step_num_1 = 1
    final_step_num_2 = 0
    
    print(f"= {final_step_num_1} - {final_step_num_2}")

    # Calculate the limit using sympy
    limit_value = sympy.limit(expression, n, sympy.oo)

    print(f"\nThe value of the limit is: {limit_value}")

solve_random_walk_problem()