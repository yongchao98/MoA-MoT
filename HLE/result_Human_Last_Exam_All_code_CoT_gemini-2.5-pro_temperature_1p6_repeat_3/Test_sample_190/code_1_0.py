import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbols for our symbolic calculation
    k, c = sympy.symbols('k c')

    # The transition probabilities for a large enough state k are defined.
    # The jump is the key value, so we map jump size to its probability.
    transitions = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
         1: sympy.Rational(1, 4) + c/k,
         2: sympy.Rational(1, 4)
    }

    print("Step 1: We define the transition probabilities based on the problem description.")
    print(f"  P(k -> k-2) = {transitions[-2]}")
    print(f"  P(k -> k-1) = {transitions[-1]}")
    print(f"  P(k -> k+1) = {transitions[1]}")
    print(f"  P(k -> k+2) = {transitions[2]}")
    print("-" * 30)

    # --- Step 2: Calculate the drift (mu_k) ---
    # The drift is the expected change in position from state k.
    # mu_k = sum(jump * probability)
    mu_k = sum(jump * prob for jump, prob in transitions.items())
    mu_k = sympy.simplify(mu_k)

    print("Step 2: Calculate the drift (mu_k), which is the expected displacement from state k.")
    mu_k_calc_str = f"(-2) * ({transitions[-2]}) + (-1) * ({transitions[-1]}) + (1) * ({transitions[1]}) + (2) * ({transitions[2]})"
    print(f"  mu_k = {mu_k_calc_str}")
    print(f"  mu_k = {mu_k}")
    print("-" * 30)

    # --- Step 3: Calculate the second moment of the jump (V_k) ---
    # This is the expected squared displacement.
    # V_k = sum(jump^2 * probability)
    V_k = sum(jump**2 * prob for jump, prob in transitions.items())
    V_k = sympy.simplify(V_k)

    print("Step 3: Calculate V_k, the expected squared displacement from state k.")
    V_k_calc_str = f"(-2)**2 * ({transitions[-2]}) + (-1)**2 * ({transitions[-1]}) + (1)**2 * ({transitions[1]}) + (2)**2 * ({transitions[2]})"
    print(f"  V_k = {V_k_calc_str}")
    print(f"  V_k = {V_k}")
    print("-" * 30)

    # --- Step 4: Apply the Lamperti/Tweedie criterion for transience ---
    # The chain is transient if lim_{k->inf} (2 * k * mu_k / V_k) > 1.
    # Let's compute this limit, which we call rho.
    # For large k, V_k approaches its limit, which is a constant.
    V_k_limit = sympy.limit(V_k, k, sympy.oo)
    rho = sympy.limit(2 * k * mu_k / V_k, k, sympy.oo)

    print("Step 4: Apply the criterion for transience.")
    print("  The chain is transient if the limit rho = lim_{k->oo} (2 * k * mu_k / V_k) is greater than 1.")
    print(f"  rho = lim_{{k->oo}} (2 * k * ({mu_k}) / ({V_k}))")
    print(f"  As k -> oo, V_k -> {V_k_limit}. So, the expression becomes:")
    print(f"  rho = lim_{{k->oo}} (2 * k * ({mu_k}) / {V_k_limit})")
    print(f"  rho = {rho}")
    print("-" * 30)

    # --- Step 5: Solve the inequality to find the condition on c ---
    # The condition for transience is rho > 1.
    # The infimum corresponds to the boundary case rho = 1.
    critical_c_solution = sympy.solve(rho - 1, c)
    critical_c = critical_c_solution[0]

    print("Step 5: Solve for c to find the transience condition.")
    print("  For the chain to be transient, we must have rho > 1:")
    print(f"  {rho} > 1")
    print(f"  Solving for c gives:")
    print(f"  c > {critical_c}")
    print("-" * 30)

    print(f"The set of values for c where the chain is transient is ({critical_c}, oo).")
    print(f"The infimum of this set is the smallest value that c can approach, which is {critical_c}.")


if __name__ == "__main__":
    solve_markov_chain_transience()