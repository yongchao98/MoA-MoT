import sympy as sp

def find_pi0():
    """
    This script symbolically derives the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define symbolic variables for the parameters
    n = sp.symbols('n', integer=True, nonnegative=True)
    rho = sp.Symbol('rho', positive=True, real=True)
    pi0 = sp.Symbol('pi_0', positive=True, real=True)

    print("Step 1: Express pi_n in terms of pi_0.")
    print("For the birth-death process with birth rate lambda_i = lambda/(i+1) and death rate mu_i = mu,")
    print("the steady-state probability pi_n can be found via the detailed balance equations.")
    print("The general solution is pi_n = pi_0 * (lambda_0*...*lambda_{n-1}) / (mu_1*...*mu_n).")
    print("Substituting the rates and defining rho = lambda/mu gives:")
    # The derived expression for pi_n is pi_0 * (rho^n / n!)
    pi_n_expr = pi0 * (rho**n / sp.factorial(n))
    print(f"pi_n = pi_0 * (rho**n / n!)\n")


    print("Step 2: Use the normalization condition Sum(pi_n for n=0 to infinity) = 1.")
    # Create the summation expression
    sum_expr = sp.Sum(pi_n_expr, (n, 0, sp.oo))
    print(f"This gives the equation: {sum_expr} = 1\n")


    print("Step 3: Evaluate the infinite sum.")
    # The sum is pi_0 * Sum(rho^n / n!), which is the Taylor series for pi_0 * e^rho
    evaluated_sum = sum_expr.doit()
    print("The summation part Sum(rho**n / n!) is the Taylor series for e^rho.")
    print(f"The equation becomes: {evaluated_sum} = 1\n")

    print("Step 4: Solve for pi_0.")
    # Create the equation object
    equation = sp.Eq(evaluated_sum, 1)
    # Solve for pi0
    solution = sp.solve(equation, pi0)

    # The solution is a list, so we take the first element
    final_pi0 = solution[0]
    
    print("The final expression for pi_0 is:")
    final_equation = sp.Eq(sp.Symbol('pi_0'), final_pi0)
    print(sp.pretty(final_equation, use_unicode=True))


find_pi0()