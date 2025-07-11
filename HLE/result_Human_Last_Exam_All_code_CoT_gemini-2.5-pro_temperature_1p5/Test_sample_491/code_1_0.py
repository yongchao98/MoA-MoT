import sympy

def solve_bargaining_problem():
    """
    Solves the bargaining problem as described.
    1. Computes the SPE payoff for Player 1.
    2. Computes the optimal payoff for Player 1 against the specified Player 2.
    3. Calculates the difference and the sum of the absolute values of the coefficients
       of the resulting rational function's numerator and denominator.
    """
    d = sympy.Symbol('d')

    # Part 1: Standard Subgame Perfect Equilibrium (SPE) Payoff for Player 1
    # Let v1 be P1's payoff, v2 be P2's payoff.
    # In equilibrium, P1 offers P2 a share (1-x) such that P2 is indifferent: 1-x = d*v2_share.
    # By symmetry, the total value to a player whose turn it is is v1_val = 1/(1+d).
    v_spe = 1 / (1 + d)
    print(f"1. Standard SPE Payoff for Player 1 (V_SPE): {v_spe}\n")

    # Part 2: Player 1's Payoff Against the Specific Player 2
    # First, calculate the expected payoff in a fully random simulation.
    # E_make: expected payoff when it's your turn to make an offer.
    # E_receive: expected payoff when it's your turn to receive an offer.
    # Offer 's' is player's own share, from U([0,1]), so E[s]=1/2.
    # E_make    = (1/2)*E[s] + (1/2)*d*E_receive  = 1/4 + (d/2)*E_receive
    # E_receive = (1/2)*E[1-s] + (1/2)*d*E_make    = 1/4 + (d/2)*E_make
    E_make, E_receive = sympy.symbols('E_make E_receive')
    eq1 = sympy.Eq(E_make, sympy.S(1)/4 + d/2 * E_receive)
    eq2 = sympy.Eq(E_receive, sympy.S(1)/4 + d/2 * E_make)
    sim_solution = sympy.solve([eq1, eq2], (E_make, E_receive))
    E_sim = sim_solution[E_make] # Payoff for the player whose turn it is to make an offer
    
    # Player 2 decides by comparing their payoff from Accepting (1-x) with their
    # expected payoff from the Reject simulation (d * E_sim).
    # P2 Accepts if 1-x >= d * E_sim
    # P1, to maximize their own share 'x', will offer the minimum to P2, so 1-x = d * E_sim.
    # P1's optimal payoff is therefore x = 1 - d * E_sim.
    v_p1 = 1 - d * E_sim
    
    print(f"2. Analyzing the specific opponent:")
    print(f"   - Expected payoff in random simulation (E_sim): {E_sim}")
    print(f"   - P1's optimal payoff against this opponent (V_P1): {sympy.simplify(v_p1)}\n")

    # Part 3: Calculate the Difference and Sum of Coefficients
    # The difference is D = V_P1 - V_SPE
    difference_expr = sympy.cancel(v_p1 - v_spe)
    
    # Get the numerator and denominator polynomials
    p_d, q_d = sympy.fraction(difference_expr)
    
    # Ensure standard form with positive leading coefficient in denominator if possible
    if q_d.as_poly(d).LC() < 0:
        p_d = -p_d
        q_d = -q_d

    print(f"3. Difference D = V_P1 - V_SPE = p(d) / q(d)")
    print(f"   The final equation is ({p_d}) / ({q_d})\n")

    p_poly = sympy.Poly(p_d, d)
    q_poly = sympy.Poly(q_d, d)

    # Get all coefficients (including zeros for missing terms)
    p_coeffs = p_poly.all_coeffs()
    q_coeffs = q_poly.all_coeffs()

    print(f"4. Calculating the sum of absolute values of coefficients:")
    print(f"   Numerator p(d) = {p_d}")
    print(f"   The integer coefficients in the numerator are: {[int(c) for c in p_coeffs]}")
    p_sum = sum(abs(c) for c in p_coeffs)
    print(f"   Sum of absolute values of numerator coefficients: {int(p_sum)}\n")
    
    print(f"   Denominator q(d) = {q_d}")
    print(f"   The integer coefficients in the denominator are: {[int(c) for c in q_coeffs]}")
    q_sum = sum(abs(c) for c in q_coeffs)
    print(f"   Sum of absolute values of denominator coefficients: {int(q_sum)}\n")
    
    total_sum = p_sum + q_sum
    print(f"5. Total sum of absolute values of coefficients = {int(p_sum)} + {int(q_sum)} = {int(total_sum)}")

if __name__ == '__main__':
    solve_bargaining_problem()