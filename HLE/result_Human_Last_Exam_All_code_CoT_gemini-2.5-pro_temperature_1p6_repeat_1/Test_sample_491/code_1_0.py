import sympy as sp

def solve_bargaining_problem():
    """
    Solves the modified Rubinstein-Stahl bargaining problem.
    """
    d = sp.Symbol('d')

    # Step 1: Calculate the payoff for Player 1 in the standard SPE
    V1_spe = 1 / (1 + d)

    # Step 2 & 3: Model the simulation world to find P2's decision threshold
    # Let e1 be P2's expected payoff (pie=1) when P1 proposes in a sim.
    # Let e2 be P2's expected payoff (pie=1) when P2 proposes in a sim.
    # In a sim, offer is U[0,1], so expected share is 0.5. A/R is 50/50.
    # e1 = 0.5 * (P2_share_if_acc) + 0.5 * (P2_payoff_if_rej)
    # e1 = 0.5 * 0.5 + 0.5 * (d * e2)
    # e2 = 0.5 * (P2_share_if_acc) + 0.5 * (P2_payoff_if_rej)
    # e2 = 0.5 * 0.5 + 0.5 * (d * e1)
    e1, e2 = sp.symbols('e1 e2')
    eq1 = sp.Eq(e1, 0.25 + 0.5 * d * e2)
    eq2 = sp.Eq(e2, 0.25 + 0.5 * d * e1)
    
    # Solve the system for e2, which is needed for P2's rejection simulation
    sol = sp.solve([eq1, eq2], (e1, e2))
    e2_val = sol[e2]
    
    # Step 4: Determine Player 1's optimal payoff against this specific P2
    # P2's threshold for accepting is based on the expected payoff from the 'Reject' simulation.
    # The 'Reject' sim starts with P2 proposing over a pie of size d.
    # Expected payoff = d * e2_val
    p2_rejection_payoff = d * e2_val
    
    # P1 will offer P2 a share of p2_rejection_payoff, keeping 1 - p2_rejection_payoff.
    # This is P1's payoff if the offer is accepted.
    v1_offer_accepted = 1 - p2_rejection_payoff
    
    # The continuation value if P1's offer is rejected is d^2 * V1_optimal.
    # P1 chooses the max of making an acceptable offer vs an unacceptable one.
    # Since d < 1, v1_offer_accepted > d^2 * v1_offer_accepted.
    # Thus, P1's optimal strategy is to make the just-acceptable offer.
    V1_optimal = v1_offer_accepted

    # Step 5: Calculate the difference and sum of absolute coefficients
    difference = V1_optimal - V1_spe

    # Simplify and extract numerator and denominator as polynomials
    difference_simplified = sp.cancel(difference)
    p_d, q_d = sp.fraction(difference_simplified)

    # Make sure coefficients are integers. Multiply by a constant if necessary.
    p_poly = sp.Poly(p_d, d)
    q_poly = sp.Poly(q_d, d)
    
    # The default sympy simplification might give rational coefficients.
    # We find the least common multiple of the denominators of the coefficients.
    p_coeffs_rat = p_poly.all_coeffs()
    q_coeffs_rat = q_poly.all_coeffs()
    
    lcm_p = sp.lcm([c.q for c in p_coeffs_rat])
    lcm_q = sp.lcm([c.q for c in q_coeffs_rat])
    common_mult = sp.lcm(lcm_p, lcm_q)
    
    p_d_int = sp.expand(p_d * common_mult)
    q_d_int = sp.expand(q_d * common_mult)

    p_poly_int = sp.Poly(p_d_int, d)
    q_poly_int = sp.Poly(q_d_int, d)
    
    p_coeffs = p_poly_int.coeffs()
    q_coeffs = q_poly_int.coeffs()

    sum_abs_coeffs_p = sum(abs(c) for c in p_coeffs)
    sum_abs_coeffs_q = sum(abs(c) for c in q_coeffs)
    total_sum_abs_coeffs = sum_abs_coeffs_p + sum_abs_coeffs_q
    
    print("Player 1's payoff in SPE:")
    print(f"V1_spe = {V1_spe}\n")
    print("Player 1's optimal payoff against the specific Player 2:")
    print(f"V1_optimal = {sp.cancel(V1_optimal)}\n")
    print("Difference (V1_optimal - V1_spe) as a quotient of polynomials p(d)/q(d):")
    print(f"p(d) = {p_d_int}")
    print(f"q(d) = {q_d_int}\n")
    print(f"Numerator p(d) coefficients are: {p_coeffs}")
    print(f"Sum of absolute values of coefficients in p(d): {sum_abs_coeffs_p}\n")
    print(f"Denominator q(d) coefficients are: {q_coeffs}")
    print(f"Sum of absolute values of coefficients in q(d): {sum_abs_coeffs_q}\n")
    print(f"Total sum of absolute values of coefficients: {total_sum_abs_coeffs}")

solve_bargaining_problem()