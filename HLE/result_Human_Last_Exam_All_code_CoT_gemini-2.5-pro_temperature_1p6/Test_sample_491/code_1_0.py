import sympy

def solve_bargaining_problem():
    """
    Solves the described bargaining game problem using symbolic mathematics.
    
    1.  Calculates the SPE payoff for Player 1.
    2.  Calculates the optimal payoff for Player 1 against the specific Player 2.
    3.  Finds the difference between these payoffs as a rational function p(d)/q(d).
    4.  Calculates the sum of the absolute values of the coefficients of p(d) and q(d).
    """
    
    # Define the discount factor 'd' as a symbolic variable
    d = sympy.symbols('d')

    # Part 1: Standard Subgame Perfect Equilibrium (SPE) Payoff for Player 1
    # In the SPE of the standard Rubinstein bargaining game, Player 1's payoff is 1 / (1 + d).
    V1_spe = 1 / (1 + d)

    # Part 2: Player 1's Optimal Payoff against the specific opponent
    # Player 2's decision to accept/reject is based on a simulation.
    # Let k be the expected fraction of the pie P2 gets in a random simulation.
    # k*pie = 0.5 * (pie/2) + 0.5 * k*d*pie
    # k = 1/4 + k*d/2  => k(1 - d/2) = 1/4 => k = 1 / (4 - 2*d)
    k = 1 / (4 - 2*d)
    
    # P2's reservation payoff is the simulated payoff from rejecting, where the pie will be 'd'.
    # Reservation_payoff_P2 = k * d
    reservation_p2 = k * d
    
    # Player 1 offers Player 2 exactly this amount to make P2 accept.
    # P1's payoff is 1 - P2's share.
    V1_optimal = 1 - reservation_p2
    
    # Part 3: Calculate the difference and the sum of coefficients
    
    # Calculate the difference between the optimal payoff and the SPE payoff
    difference_expr = V1_optimal - V1_spe
    
    # Simplify the expression to get a single fraction p(d)/q(d)
    # sympy.cancel() is used for rational functions to put them into p/q form.
    final_fraction = sympy.cancel(difference_expr)
    
    p_d, q_d = sympy.fraction(final_fraction)
    
    # Ensure numerator and denominator are expanded polynomials
    p_d_expanded = sympy.expand(p_d)
    q_d_expanded = sympy.expand(q_d)

    # Represent as polynomial objects to extract coefficients
    p_poly = sympy.Poly(p_d_expanded, d)
    q_poly = sympy.Poly(q_d_expanded, d)
    
    p_coeffs = p_poly.all_coeffs()
    q_coeffs = q_poly.all_coeffs()

    # Calculate the sum of the absolute values of the coefficients
    sum_abs_coeffs_p = sum(abs(c) for c in p_coeffs)
    sum_abs_coeffs_q = sum(abs(c) for c in q_coeffs)
    total_sum = sum_abs_coeffs_p + sum_abs_coeffs_q
    
    print("The difference between the two payoffs is p(d)/q(d), where:")
    print(f"p(d) = {p_d_expanded}")
    print(f"q(d) = {q_d_expanded}")
    print("\nThe coefficients of the numerator p(d) are:", [int(c) for c in p_coeffs])
    print("The coefficients of the denominator q(d) are:", [int(c) for c in q_coeffs])
    print("\nThe sum of the absolute values of these coefficients is:")
    print(f"(|{p_coeffs[0]}| + |{p_coeffs[1]}|) + (|{q_coeffs[0]}| + |{q_coeffs[1]}| + |{q_coeffs[2]}|) = {int(sum_abs_coeffs_p)} + {int(sum_abs_coeffs_q)} = {int(total_sum)}")

solve_bargaining_problem()