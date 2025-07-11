import sympy

def solve_bargaining_problem():
    """
    Solves the bargaining problem as described.
    1. Computes the SPE payoff for Player 1.
    2. Computes the optimal payoff for Player 1 against the specific Player 2.
    3. Calculates the difference and the sum of absolute values of coefficients.
    """
    # Define d as a symbol
    d = sympy.Symbol('d')

    # Part 1: Standard Subgame Perfect Equilibrium (SPE) Payoff
    # v1_spe is P1's payoff when P1 proposes
    # v2_spe is P1's payoff when P2 proposes
    v1_spe = sympy.Function('v1_spe')(d)
    v2_spe = sympy.Function('v2_spe')(d)

    # System of equations for SPE:
    # eq1: v1 = 1 - d*(1-v2) where 1-v2 is P2's payoff when P2 proposes
    # eq2: v2 = d*v1 where P1 accepts the minimum offer from P2
    eq1_spe = sympy.Eq(v1_spe, 1 - d + d * v2_spe)
    eq2_spe = sympy.Eq(v2_spe, d * v1_spe)

    # Solve for v1_spe
    sol_spe = sympy.solve([eq1_spe, eq2_spe], [v1_spe, v2_spe])
    V1_SPE = sympy.cancel(sol_spe[v1_spe])

    # Part 2: Payoff against the specific Player 2
    # E1 is P2's expected payoff when P1 (sim) proposes
    # E2 is P2's expected payoff when P2 (sim) proposes
    E1 = sympy.Function('E1')(d)
    E2 = sympy.Function('E2')(d)

    # System of equations for the random simulation expected payoffs.
    # A random offer z~U[0,1] has E[z]=0.5. A random A/R is 50/50.
    # E1 = 0.5 * E[1-z] + 0.5 * d*E2 = 0.5*0.5 + 0.5*d*E2 = 1/4 + d*E2/2
    # E2 = 0.5 * E[1-z] + 0.5 * d*E1 = 0.5*0.5 + 0.5*d*E1 = 1/4 + d*E1/2
    eq_E1 = sympy.Eq(E1, sympy.Rational(1, 4) + sympy.Rational(1, 2) * d * E2)
    eq_E2 = sympy.Eq(E2, sympy.Rational(1, 4) + sympy.Rational(1, 2) * d * E1)

    # Solve for E2
    sol_E = sympy.solve([eq_E1, eq_E2], [E1, E2])
    E2_val = sympy.cancel(sol_E[E2])

    # P2 accepts P1's offer x if 1-x is greater than P2's expected payoff from rejecting.
    # P2's expected payoff from rejecting is d * E2_val.
    # P1 plays optimally by offering x such that 1-x = d * E2_val.
    # P1's payoff is therefore x.
    V1_optimal = sympy.cancel(1 - d * E2_val)

    # Part 3: Calculate the difference and the sum of coefficients
    Difference = sympy.cancel(V1_optimal - V1_SPE)

    # Extract numerator and denominator polynomials
    p_d, q_d = sympy.fraction(Difference)

    p_poly = sympy.Poly(p_d, d)
    q_poly = sympy.Poly(q_d, d)

    # For a canonical form, ensure the leading coefficient of the denominator is positive
    if q_poly.LC() < 0:
        p_poly = -p_poly
        q_poly = -q_poly
    
    p_coeffs = [int(c) for c in p_poly.all_coeffs()]
    q_coeffs = [int(c) for c in q_poly.all_coeffs()]

    # Calculate sum of absolute values
    p_sum_abs = sum(abs(c) for c in p_coeffs)
    q_sum_abs = sum(abs(c) for c in q_coeffs)
    total_sum = p_sum_abs + q_sum_abs

    # --- Print the results ---
    print(f"Player 1's payoff in the standard SPE game: V1_SPE = {V1_SPE}")
    print(f"Player 1's optimal payoff against the specific opponent: V1_optimal = {V1_optimal}")
    print("-" * 50)
    print(f"The difference is p(d)/q(d) = {p_poly.as_expr()} / ({q_poly.as_expr()})")
    print("-" * 50)

    # Print numerator details
    print(f"Numerator polynomial is p(d) = {p_poly.as_expr()}")
    p_terms = []
    for i, c in enumerate(p_coeffs):
        power = p_poly.degree() - i
        if c != 0:
            p_terms.append(f"({c})*d**{power}")
    p_eq_str = " + ".join(p_terms)
    print(f"Numerator equation: p(d) = {p_eq_str}")
    print(f"The integer coefficients are: {p_coeffs}")
    print(f"Sum of absolute values of numerator coefficients: {p_sum_abs}")
    print("-" * 50)
    
    # Print denominator details
    print(f"Denominator polynomial is q(d) = {q_poly.as_expr()}")
    q_terms = []
    for i, c in enumerate(q_coeffs):
        power = q_poly.degree() - i
        # Handle constant term power display
        if power > 1:
          q_terms.append(f"({c})*d**{power}")
        elif power == 1:
          q_terms.append(f"({c})*d")
        else:
          q_terms.append(f"({c})")

    q_eq_str = " + ".join(q_terms)
    print(f"Denominator equation: q(d) = {q_eq_str}")
    print(f"The integer coefficients are: {q_coeffs}")
    print(f"Sum of absolute values of denominator coefficients: {q_sum_abs}")
    print("-" * 50)

    print(f"Total sum of absolute values of coefficients = {p_sum_abs} + {q_sum_abs} = {total_sum}")

if __name__ == '__main__':
    solve_bargaining_problem()
