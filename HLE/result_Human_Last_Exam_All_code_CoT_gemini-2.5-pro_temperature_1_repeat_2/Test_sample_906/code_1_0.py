import sympy
from sympy import Sum, factorial, oo, exp, Eq, solve

def solve_steady_state_probability():
    """
    This function symbolically derives the steady-state probability pi_0
    for a birth-death process with rates lambda_i = lambda / (i + 1) and mu_i = mu.
    """

    # Define the symbolic variables used in the derivation.
    k = sympy.Symbol('k', integer=True, nonnegative=True)
    rho = sympy.Symbol('rho', real=True, positive=True)
    pi_0_sym = sympy.Symbol('pi_0', real=True, positive=True)

    print("Step 1: Express the steady-state probability pi_k in terms of pi_0 and rho.")
    print("From the detailed balance equations, we find the relation:")
    # The general expression for pi_k in terms of pi_0 is pi_0 * (rho^k / k!).
    pi_k_expression = pi_0_sym * rho**k / factorial(k)
    print(f"pi_k = {pi_k_expression}")
    print("-" * 50)

    print("Step 2: Use the normalization condition Sum(pi_k for k=0 to inf) = 1.")
    # The sum of all probabilities must equal 1.
    total_probability_sum = Sum(pi_k_expression, (k, 0, oo))
    normalization_equation = Eq(total_probability_sum, 1)
    print("The normalization equation is:")
    print(normalization_equation)
    print("-" * 50)

    print("Step 3: Evaluate the infinite sum.")
    # The sum is the Taylor series for exp(rho). Sympy can evaluate it directly.
    evaluated_sum = total_probability_sum.doit()
    evaluated_equation = Eq(evaluated_sum, 1)
    print("The sum Sum(rho**k / k!, (k, 0, oo)) is the series for exp(rho).")
    print("The equation simplifies to:")
    print(evaluated_equation)
    print("-" * 50)

    print("Step 4: Solve the simplified equation for pi_0.")
    # Solve the equation pi_0 * exp(rho) = 1 for pi_0.
    solution = solve(evaluated_equation, pi_0_sym)
    final_expression = solution[0]
    
    # Present the final formula clearly.
    final_equation = Eq(sympy.Symbol('pi_0'), final_expression)
    print("The final expression for pi_0 is:")
    print(f"{final_equation.lhs} = {final_equation.rhs}")
    print("\nIn this final equation:")
    print(f"* The term 'pi_0' is the steady-state probability of being in state 0.")
    print(f"* The term 'rho' is the ratio of the base arrival rate to the service rate (lambda / mu).")
    print(f"* The term 'exp' is the exponential function (e).")

if __name__ == '__main__':
    solve_steady_state_probability()