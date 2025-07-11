import sympy as sp

def solve_steady_state():
    """
    Symbolically derives the steady-state probability pi_0 for the given
    birth-death process.
    """
    # Define symbols for the parameters
    rho = sp.symbols('rho', real=True, positive=True)
    pi_0_sym = sp.symbols('pi_0')
    k = sp.symbols('k', integer=True, nonnegative=True)

    # From the detailed balance equations, we derive the general form for pi_k
    # pi_k = (rho^k / k!) * pi_0
    pi_k_relation = (rho**k / sp.factorial(k)) * pi_0_sym

    print("Step 1: The relationship between pi_k and pi_0 is derived from the balance equations.")
    print(f"pi_k = (rho^k / k!) * pi_0")
    print("-" * 30)

    # Step 2: Use the normalization condition Sum(pi_k) = 1 for k from 0 to infinity
    # This leads to pi_0 * Sum(rho^k / k!) = 1
    sum_term = pi_k_relation / pi_0_sym
    
    # The sum is the Taylor series for exp(rho)
    series_sum = sp.Sum(sum_term, (k, 0, sp.oo)).doit()

    print("Step 2: The normalization condition is Sum(pi_k) for k=0 to inf = 1.")
    print("This means: pi_0 * Sum(rho^k / k!) = 1")
    print(f"The sum Sum(rho^k / k!) evaluates to: {series_sum}")
    print("-" * 30)

    # Step 3: Solve for pi_0
    # We set up the equation: pi_0 * exp(rho) = 1
    normalization_eq = sp.Eq(pi_0_sym * series_sum, 1)
    
    # Solve the equation for pi_0
    pi_0_solution = sp.solve(normalization_eq, pi_0_sym)

    # The result is a list containing the solution
    final_expression = pi_0_solution[0]

    print("Step 3: Solve the equation for pi_0.")
    print(f"Equation: {pi_0_sym} * {series_sum} = 1")
    
    # We print each part of the final equation for clarity.
    print("\nFinal equation:")
    print(f"pi_0 = {sp.pretty(final_expression)}")


solve_steady_state()