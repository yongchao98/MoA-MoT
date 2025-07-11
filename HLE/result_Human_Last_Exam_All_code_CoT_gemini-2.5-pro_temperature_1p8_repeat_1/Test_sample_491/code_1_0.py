import sympy as sp

def solve_bargaining_problem():
    """
    Solves the bargaining problem as described.
    1. Calculates the SPE payoff for Player 1.
    2. Calculates Player 1's payoff against the specified Player 2.
    3. Finds the difference and expresses it as p(d)/q(d).
    4. Calculates the sum of the absolute values of the coefficients of p(d) and q(d).
    """
    d = sp.Symbol('d')

    # Step 1: Payoff in the unique subgame perfect equilibrium (SPE)
    # Let V1 be P1's payoff when P1 offers, V2 be P1's payoff when P2 offers.
    # P1 offers x. P2 accepts if 1-x >= d*(1-V2'), where V2' is P2's payoff in next round.
    # In SPE, P2's payoff V2' equals 1-V1 where V1 is the value to P1.
    # This leads to P1 offering x s.t. 1-x = d*(1-(1/(1+d))) * (size of pie) oh wait.
    # A simpler formulation: V1 = 1 - d*V_P2_offering. P2's total value is (1-V2).
    # Correct SPE derivation:
    # V1 = 1 - d * (value to P2 if P2 offers). P2 gets 1-V2.
    # But V2 (P1's payoff) = d*V1.
    # So P2 gets 1-d*V1.
    # And value to P2 if P2 offers is 1-V2 = 1-d*V1. Wait, this isn't right.
    # The value to P2 in the subgame where P2 offers is 1-V2 (since V2 is P1's value).
    # P1 offers x, P2 gets 1-x. P2 rejects if d*(1-V2) > 1-x. P1 offers min x s.t. 1-x = d*(1-V2).
    # V1 = x = 1 - d*(1-V2).
    # V2 (P1 payoff when P2 offers) = d*V1.
    # V1 = 1 - d*(1 - d*V1) = 1 - d + d**2*V1 => V1(1-d**2) = 1-d => V1 = (1-d)/((1-d)(1+d))
    P_spe = 1 / (1 + d)

    # Step 2: Payoff against the specified Player 2
    # P2's A/R decision is based on a simulation.
    # Let E_s1 be P2's expected share of the pie in sim when P1 offers.
    # Let E_s2 be P2's expected share of the pie in sim when P2 offers.
    # In simulation, choices are random (1/2 prob) and offers are U[0,1].
    # Expected share P2 gets from an offer is E[1-s] = 0.5.
    E_s1, E_s2 = sp.symbols('E_s1 E_s2')
    eq1 = sp.Eq(E_s1, 0.5 * 0.5 + 0.5 * d * E_s2)  # 0.5*accept(get 0.5) + 0.5*reject(get d*E_s2)
    eq2 = sp.Eq(E_s2, 0.5 * 0.5 + 0.5 * d * E_s1)  # 0.5*accept(get 0.5) + 0.5*reject(get d*E_s1)
    
    solution = sp.solve([eq1, eq2], (E_s1, E_s2))
    simulated_share_p2 = solution[E_s2]

    # When P1 makes an offer, P2 compares accepting with the simulated payoff of rejecting.
    # Rejection leads to a subgame where pie is d and it's P2's turn.
    # P2's simulated rejection payoff = d * simulated_share_p2
    rejection_payoff = d * simulated_share_p2
    
    # P1 wants to offer x to maximize x, subject to P2 accepting.
    # P2 accepts if P2's share (1-x) >= rejection_payoff.
    # So P1 sets 1-x = rejection_payoff => x = 1 - rejection_payoff.
    P_mod = 1 - rejection_payoff

    # Step 3: Compute the difference
    difference = P_mod - P_spe

    # Step 4: Express as p(d)/q(d) and find coefficients
    # Use sp.together to combine into a single fraction
    difference_fraction = sp.together(difference)
    
    # Extract numerator and denominator and expand them to standard polynomial form
    p_d_raw, q_d_raw = difference_fraction.as_numer_denom()
    p_d = sp.expand(p_d_raw)
    q_d = sp.expand(q_d_raw)
    
    # Extract coefficients from the polynomials. all_coeffs() includes zero terms.
    p_coeffs = sp.Poly(p_d, d).all_coeffs()
    q_coeffs = sp.Poly(q_d, d).all_coeffs()
    
    # Calculate the sum of the absolute values of the coefficients
    sum_abs_coeffs = sum(abs(c) for c in p_coeffs) + sum(abs(c) for c in q_coeffs)
    
    # Step 5: Print the results
    print("The calculated difference, p(d)/q(d), is:")
    # Using sp.printing.python() to get a parseable string expression
    # The requested output wants "each number in the final equation"
    # So we format the output to be clear.
    p_d_str = str(p_d).replace('**', '^')
    q_d_str = str(q_d).replace('**', '^')
    print(f"\nFinal Equation:\n")
    print(f"Difference = ({p_d_str}) / ({q_d_str})")

    # Outputting each coefficient ("number") in the equation explicitly
    print(f"\nNumerator p(d) = {p_d_str}")
    print(f"Numerator Coefficients: {p_coeffs}")
    print(f"Denominator q(d) = {q_d_str}")
    print(f"Denominator Coefficients: {q_coeffs}")

    print(f"\nThe sum of the absolute values of the integer coefficients is:")
    print(int(sum_abs_coeffs))


solve_bargaining_problem()