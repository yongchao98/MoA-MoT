import sympy

def solve_probability():
    """
    Calculates the upper bound for the probability that the production process
    reaches a point where exactly 50% of the products are good.
    """
    # Define the symbolic variable for our probability distribution
    x = sympy.Symbol('x')

    # Initial number of good (W_0) and defective (B_0) products
    W0 = 2
    B0 = 1

    # The limit of the proportion of good products, X_infinity, follows a Beta distribution
    # Beta(W_0, B_0). The Probability Density Function (PDF) for Beta(2, 1) is f(x) = 2x.
    pdf = 2 * x

    # The process stops when the proportion of good products is 1/2.
    # If the process never stops (T=infinity), the proportion must always be > 1/2.
    # Thus, we need to calculate the conditional expectation E[X_infinity | X_infinity >= 1/2].

    # 1. Calculate P(X >= 1/2), the probability of the condition
    prob_ge_half = sympy.integrate(pdf, (x, sympy.S.Half, 1))

    # 2. Calculate the numerator for the conditional expectation: E[X * 1_{X >= 1/2}]
    integrand = x * pdf
    exp_x_numerator = sympy.integrate(integrand, (x, sympy.S.Half, 1))

    # 3. Calculate the conditional expectation E[X_infinity | T = infinity]
    E_X_inf_if_T_inf = exp_x_numerator / prob_ge_half

    # Now, we use the Optional Stopping Theorem for the martingale X_t.
    # Let p = P(T < infinity).
    # E[X_0] = p * E[X_T | T < infinity] + (1-p) * E[X_infinity | T = infinity]

    # Initial expectation E[X_0]
    E_X0 = sympy.S(W0) / (W0 + B0)

    # Expectation at stopping time T, E[X_T | T < infinity]
    E_X_T_if_T_finite = sympy.S.Half

    # We have all the components for the equation. Now, we solve for p.
    # The equation is: E_X0 = p * E_X_T_if_T_finite + (1-p) * E_X_inf_if_T_inf
    # Rearranging to solve for p gives:
    # p = (E_X_inf_if_T_inf - E_X0) / (E_X_inf_if_T_inf - E_X_T_if_T_finite)
    p = (E_X_inf_if_T_inf - E_X0) / (E_X_inf_if_T_inf - E_X_T_if_T_finite)

    # The calculated probability 'p' is the least upper bound.
    # Print the steps of the calculation and the final equation.
    print("--- Calculation Steps ---")
    print(f"Initial state: W_0 = {W0}, B_0 = {B0}")
    print(f"Martingale X_t = W_t / (W_t + B_t)")
    print(f"Initial value E[X_0] = {W0}/{W0+B0} = {E_X0}")
    print(f"Stopping value X_T = 1/2 = {E_X_T_if_T_finite}")
    print("\nAnalysis for non-stopping case (T=infinity):")
    print("X_infinity ~ Beta(2, 1) with PDF f(x) = 2x, and X_infinity >= 1/2")
    print(f"E[X_infinity | T=infinity] = {E_X_inf_if_T_inf}")
    print("\n--- Final Equation (from Optional Stopping Theorem) ---")
    print("Let p = P(T < infinity). The equation is:")
    print(f"E[X_0] = p * E[X_T | T < inf] + (1-p) * E[X_inf | T = inf]")
    print(f"This gives the concrete equation:")
    print(f"{E_X0} = p * {E_X_T_if_T_finite} + (1-p) * {E_X_inf_if_T_inf}")
    print("\n--- Result ---")
    print(f"The calculated probability, which is the least upper bound, is p = {p}")

solve_probability()