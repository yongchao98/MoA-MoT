import sympy

def solve_bargaining_problem():
    """
    This function solves the modified Rubinstein-Stahl bargaining problem.
    It follows the plan outlined above to find the required sum of absolute coefficients.
    """

    # Define the discount factor 'd' as a symbolic variable
    d = sympy.Symbol('d')

    # Step 1: Standard Subgame Perfect Equilibrium (SPE) Payoff for Player 1
    # In a standard game, Player 1's payoff is s1 = 1 / (1 + d).
    payoff_spe = 1 / (1 + d)

    # Step 2: Calculate Player 2's simulated rejection value.
    # The simulated value is the expected share for the offering player in a fully random game.
    # Let v_o be this share. The recursive equation is v_o = 1/4 + (d/2) * (1/4 + (d/2)*v_o),
    # which solves to v_o = 1 / (2 * (2 - d)).
    simulated_share_if_offer = 1 / (2 * (2 - d))
    
    # If P2 rejects P1's initial offer, the pie size becomes d and it's P2's turn.
    # So the value P2 uses for comparison is d * the simulated share.
    p2_rejection_threshold = d * simulated_share_if_offer

    # Step 3: Determine Player 1's optimal payoff in the modified game.
    # P1 offers P2 exactly the threshold value to secure acceptance, so P1 gets 1 minus that value.
    # We verified that P1 making an unacceptable offer is a worse strategy.
    payoff_optimal = 1 - p2_rejection_threshold

    # Step 4: Calculate the difference between the optimal and SPE payoffs.
    difference = payoff_optimal - payoff_spe

    # Step 5: Simplify the expression to the form p(d)/q(d).
    difference_simplified = sympy.cancel(difference)
    p_d, q_d = sympy.fraction(difference_simplified)
    
    # Re-order polynomials to standard form (highest power first).
    p_d = sympy.Poly(p_d, d).as_expr()
    q_d = sympy.Poly(q_d, d).as_expr()
    
    # Step 6: Extract coefficients and compute the sum of their absolute values.
    p_coeffs = sympy.Poly(p_d, d).all_coeffs()
    q_coeffs = sympy.Poly(q_d, d).all_coeffs()

    sum_abs_coeffs = sum(abs(c) for c in p_coeffs) + sum(abs(c) for c in q_coeffs)

    # --- Output the results as requested ---
    
    # The final equation p(d)/q(d)
    print(f"The calculated difference p(d)/q(d) is:\n")
    sympy.init_printing(use_unicode=False)
    final_eq_str = f"({sympy.printing.sstr(p_d)}) / ({sympy.printing.sstr(q_d)})"
    print(final_eq_str)
    
    # The numbers in the final equation (the coefficients)
    p_final_coeffs_int = [int(c) for c in p_coeffs]
    q_final_coeffs_int = [int(c) for c in q_coeffs]
    
    print("\n--- Equation Numbers Breakdown ---")
    print(f"Numbers (coefficients) in the numerator p(d): {p_final_coeffs_int}")
    print(f"Numbers (coefficients) in the denominator q(d): {q_final_coeffs_int}")

    # The final sum
    print("\n--- Final Answer ---")
    print(f"The sum of the absolute values of these coefficients is: {int(sum_abs_coeffs)}")

solve_bargaining_problem()