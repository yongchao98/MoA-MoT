import sympy as sp

def solve_steady_state_pi0():
    """
    This function symbolically derives the steady-state probability pi_0
    for the given birth-death process.
    """
    # Step 1: Define symbolic variables
    # rho represents lambda/mu
    # i is the state index
    # pi_0 is the steady-state probability of being in state 0
    rho = sp.symbols('rho', positive=True)
    i = sp.symbols('i', integer=True, nonnegative=True)
    pi_0 = sp.symbols('pi_0')

    print("Step 1: Express pi_i in terms of pi_0 and rho.")
    print("From the detailed balance equations for the given rates, we derive the general formula:")
    # Using sp.pprint for a nicer mathematical output
    pi_i_expr = pi_0 * (rho**i / sp.factorial(i))
    print(f"pi_i = {sp.printing.pretty(pi_i_expr, use_unicode=False)}")
    print("-" * 30)

    print("Step 2: Use the normalization condition Sum(pi_i) = 1 for i from 0 to infinity.")
    print("This gives the equation: pi_0 * Sum(rho**i / i!) = 1")
    print("-" * 30)

    # Step 3: Evaluate the infinite sum
    # The term inside the sum is rho**i / i!
    series_term = rho**i / sp.factorial(i)
    # The sum is the Taylor series for e^rho
    infinite_sum = sp.summation(series_term, (i, 0, sp.oo))

    print("Step 3: Evaluate the infinite sum Sum(rho**i / i!).")
    print(f"The sum is the Taylor series for e^rho, which SymPy evaluates to: {infinite_sum}")
    print("-" * 30)

    # Step 4: Solve for pi_0
    # The equation becomes pi_0 * infinite_sum = 1
    equation = sp.Eq(pi_0 * infinite_sum, 1)
    solution = sp.solve(equation, pi_0)
    final_expression = solution[0]

    print("Step 4: Solve the equation for pi_0.")
    print(f"The equation is: {sp.printing.pretty(equation, use_unicode=False)}")
    print("\nThe final expression for pi_0 is:")
    # We want to show the final equation clearly
    final_eq_str = f"pi_0 = 1 / {infinite_sum}"
    print(final_eq_str)
    print(f"pi_0 = {final_expression}")


if __name__ == '__main__':
    solve_steady_state_pi0()