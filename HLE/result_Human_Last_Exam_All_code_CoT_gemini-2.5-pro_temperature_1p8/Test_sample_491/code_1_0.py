import sympy

def solve_bargaining_problem():
    """
    Solves the specified bargaining problem by performing the following steps:
    1.  Calculates Player 2's decision threshold based on simulation.
    2.  Determines Player 1's optimal payoff against this specific Player 2.
    3.  Calculates Player 1's payoff in the standard SPE.
    4.  Finds the difference and expresses it as p(d)/q(d).
    5.  Calculates the sum of the absolute values of the coefficients of p(d) and q(d).
    """
    # Define d as a symbol
    d = sympy.Symbol('d')

    # Step 1 & 2: Calculate Player 2's expected payoff in the 'Reject' simulation
    # In the simulation, when P2 rejects P1's initial offer, the pie size is d
    # and it becomes P2's turn. We need to find P2's expected payoff from this point.
    # Let E2_i(S) be Player 2's expected simulated payoff when the pie is size S
    # and it is Player i's turn to offer. We are looking for E2_2(d).
    #
    # The system of equations for the expected payoff coefficients (k_i = E2_i(S)/S) is:
    # k2 = 1/4 + (d/2) * k1  (When it's P2's turn)
    # k1 = 1/4 + (d/2) * k2  (When it's P1's turn)
    k1, k2 = sympy.symbols('k1 k2')
    eq1 = sympy.Eq(k2, sympy.Rational(1, 4) + d/2 * k1)
    eq2 = sympy.Eq(k1, sympy.Rational(1, 4) + d/2 * k2)
    sol = sympy.solve([eq1, eq2], (k1, k2))
    k2_sol = sol[k2]  # This is E2_2(S)/S

    # The expected payoff from the 'Reject' simulation is E2_2(d) = d * k2_sol
    payoff_R_sim = d * k2_sol
    
    # This threshold is the minimum offer x* Player 1 must make.
    x_star = sympy.simplify(payoff_R_sim)

    # Step 3: Player 1's optimal payoff is 1 - x*
    payoff_P1_optimal = 1 - x_star
    payoff_P1_optimal = sympy.simplify(payoff_P1_optimal)

    # Step 4: Calculate Player 1's payoff in the standard SPE
    payoff_P1_spe = 1 / (1 + d)

    # Step 5: Compute the difference
    difference = sympy.simplify(payoff_P1_optimal - payoff_P1_spe)

    # Step 6: Express as p(d)/q(d) and sum coefficients
    p_d, q_d = sympy.fraction(difference)

    # Expand polynomials to ensure all terms are present
    p_d_expanded = sympy.expand(p_d)
    q_d_expanded = sympy.expand(q_d)

    print("The difference between the two payoffs is the rational function p(d)/q(d):")
    # Output the numbers in the final equation by printing the polynomials
    print(f"p(d) = {p_d_expanded}")
    print(f"q(d) = {q_d_expanded}")

    # Extract coefficients
    poly_p = sympy.Poly(p_d_expanded, d)
    poly_q = sympy.Poly(q_d_expanded, d)
    
    coeffs_p = [int(c) for c in poly_p.all_coeffs()]
    coeffs_q = [int(c) for c in poly_q.all_coeffs()]

    print("\nThe integer coefficients for the numerator and denominator are:")
    print(f"Numerator p(d): {coeffs_p}")
    print(f"Denominator q(d): {coeffs_q}")

    # Calculate sum of absolute values
    sum_abs_coeffs_p = sum(abs(c) for c in coeffs_p)
    sum_abs_coeffs_q = sum(abs(c) for c in coeffs_q)
    total_sum = sum_abs_coeffs_p + sum_abs_coeffs_q

    print(f"\nThe sum of absolute values for the numerator is: {sum_abs_coeffs_p}")
    print(f"The sum of absolute values for the denominator is: {sum_abs_coeffs_q}")
    print(f"The total sum of the absolute values of the coefficients is: {total_sum}")

if __name__ == '__main__':
    solve_bargaining_problem()