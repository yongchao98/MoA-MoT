import sympy

def solve_bargaining_problem():
    """
    Solves the bargaining problem as described.
    1. Calculates the SPE payoff for Player 1 in the standard game.
    2. Calculates the optimal payoff for Player 1 against the specified opponent.
    3. Computes the difference and the sum of the absolute values of the coefficients
       of the resulting rational function's numerator and denominator.
    """
    d = sympy.Symbol('d')

    # Part 1: Standard Subgame Perfect Equilibrium (SPE) Payoff for Player 1
    # Let v1 and v2 be the expected payoffs for P1 and P2 when it's their turn.
    # v1 = 1 - d*v2
    # v2 = 1 - d*v1
    # Solving this system gives v1 = 1 / (1 + d)
    payoff_spe = 1 / (1 + d)

    # Part 2: Payoff against the specific Player 2
    # First, calculate Player 2's expected payoff in the 'Reject' simulation.
    # Let v_sim1 and v_sim2 be the expected payoff shares for Player 2 when it's
    # Player 1's and Player 2's turn to propose in the simulation, respectively.
    # In the simulation: offers are s~U[0,1], A/R is 50/50.
    # E[s] = E[1-s] = 1/2.
    # When P1 proposes, P2 is offered (1-s)*pie. If P2 accepts (50%), payoff is (1-s)*pie.
    # If P2 rejects (50%), payoff is v_sim2*d*pie.
    # v_sim1*pie = 0.5 * E[(1-s)*pie] + 0.5 * (v_sim2 * d * pie) => v_sim1 = 0.5*0.5 + 0.5*d*v_sim2 = 1/4 + d/2*v_sim2
    # When P2 proposes, P2 keeps s*pie. If P1 accepts (50%), payoff is s*pie.
    # If P1 rejects (50%), payoff is v_sim1*d*pie.
    # v_sim2*pie = 0.5 * E[s*pie] + 0.5 * (v_sim1 * d*pie) => v_sim2 = 0.5*0.5 + 0.5*d*v_sim1 = 1/4 + d/2*v_sim1
    v_sim1 = sympy.Symbol('v_sim1')
    v_sim2 = sympy.Symbol('v_sim2')
    
    eq1 = sympy.Eq(v_sim1, sympy.Rational(1, 4) + d/2 * v_sim2)
    eq2 = sympy.Eq(v_sim2, sympy.Rational(1, 4) + d/2 * v_sim1)
    
    # Solve the system for v_sim2
    sim_solution = sympy.solve([eq1, eq2], (v_sim1, v_sim2))
    v_sim2_val = sim_solution[v_sim2]
    
    # Player 2's expected payoff in the 'Reject' simulation (pie starts at size d)
    reject_payoff = d * v_sim2_val

    # Player 1's optimal offer x makes Player 2 indifferent.
    # P2's payoff from accepting = 1 - x
    # P2 accepts if 1 - x >= reject_payoff.
    # P1 maximizes x, so x = 1 - reject_payoff.
    payoff_optimal = 1 - reject_payoff
    
    # Part 3: Compute the difference and the sum of coefficients
    difference = sympy.cancel(payoff_optimal - payoff_spe)

    # Get numerator and denominator
    p_d, q_d = sympy.fraction(difference)
    
    # Ensure leading coefficients are positive if possible (standard form)
    # Note: sympy.fraction might return them with negative leading coeffs.
    # The sum of absolute values will be the same regardless.
    if q_d.as_poly(d).LC() < 0:
        p_d = -p_d
        q_d = -q_d

    p_poly = sympy.Poly(p_d, d)
    q_poly = sympy.Poly(q_d, d)

    p_coeffs = [int(c) for c in p_poly.all_coeffs()]
    q_coeffs = [int(c) for c in q_poly.all_coeffs()]
    
    p_coeffs_abs_sum = sum(abs(c) for c in p_coeffs)
    q_coeffs_abs_sum = sum(abs(c) for c in q_coeffs)
    
    total_sum = p_coeffs_abs_sum + q_coeffs_abs_sum
    
    print("--- Calculation Steps ---")
    print(f"Payoff 1 (Optimal vs Specific Opponent): {sympy.simplify(payoff_optimal)}")
    print(f"Payoff 2 (Standard SPE): {payoff_spe}")
    print(f"Difference p(d)/q(d) = ({p_d}) / ({q_d})")
    print("\n--- Final Equation Coefficients ---")
    print(f"Numerator p(d) has coefficients: {p_coeffs}")
    print(f"Denominator q(d) has coefficients: {q_coeffs}")
    print(f"Sum of absolute values of numerator coefficients: {p_coeffs_abs_sum}")
    print(f"Sum of absolute values of denominator coefficients: {q_coeffs_abs_sum}")
    print("\n--- Final Answer ---")
    print(f"The total sum of the absolute values of the coefficients is: {total_sum}")
    
solve_bargaining_problem()
<<<14>>>