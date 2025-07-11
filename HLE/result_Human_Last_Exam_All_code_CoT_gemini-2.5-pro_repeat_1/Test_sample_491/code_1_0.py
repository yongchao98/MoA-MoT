import sympy

def solve_bargaining_problem():
    """
    Solves the described bargaining problem by:
    1. Calculating the SPE payoff for Player 1.
    2. Calculating Player 1's optimal payoff against the specified Player 2.
    3. Finding the difference and expressing it as a ratio of polynomials p(d)/q(d).
    4. Summing the absolute values of the integer coefficients of p(d) and q(d).
    """
    d = sympy.Symbol('d')

    # Part 1: Standard Subgame Perfect Equilibrium (SPE) Payoff for Player 1
    # In the SPE, P1 offers x such that P2 is indifferent between accepting 1-x
    # and rejecting to get d*V1_SPE (P2's payoff in the next round, by symmetry).
    # So, 1-x = d*V1_SPE. P1's payoff is x, so V1_SPE = 1 - d*V1_SPE.
    # V1_SPE * (1 + d) = 1  => V1_SPE = 1 / (1 + d)
    V1_SPE = 1 / (1 + d)

    # Part 2: Payoff against the specific Player 2
    # First, calculate the expected payoff for P2 in the rejection simulation.
    # In the simulation, play is random. Let E_P1_offers and E_P2_offers be P2's
    # expected share of a pie of size 1 when P1 or P2 offers, respectively.
    E_P1_offers = sympy.Symbol('E_P1_offers')
    E_P2_offers = sympy.Symbol('E_P2_offers')

    # The expected value of a U[0,1] offer for the other player is 0.5.
    # The expected value of keeping the rest is 1 - 0.5 = 0.5.
    # With 1/2 prob of accept/reject:
    eq1 = sympy.Eq(E_P1_offers, (1/2) * (1 - 1/2) + (1/2) * d * E_P2_offers)
    eq2 = sympy.Eq(E_P2_offers, (1/2) * (1 - 1/2) + (1/2) * d * E_P1_offers)
    
    # Solve this system for E_P2_offers
    sim_payoff_solution = sympy.solve([eq1, eq2], (E_P1_offers, E_P2_offers))
    E_P2_sim_payoff = sim_payoff_solution[E_P2_offers]

    # Player 2 will reject P1's offer of x if their share (1-x) is less than
    # their simulated payoff from rejecting. The simulated payoff is for the next round,
    # where the pie is size d and P2 is the offerer.
    rejection_threshold = d * E_P2_sim_payoff
    
    # P1's optimal strategy is to offer x such that P2's share 1-x is exactly
    # at the rejection threshold. This maximizes P1's share x while ensuring acceptance.
    # We must check that this is better than making a rejectable offer.
    # Payoff from acceptable offer: V1_accept = 1 - rejection_threshold
    # Payoff from rejectable offer: P2 offers 0, P1 rejects, game continues to t=2.
    # Payoff is d^2 * V1_special. Since d < 1, d^2 * V < V. The acceptable offer is optimal.
    V1_special = 1 - rejection_threshold
    
    # Part 3: Compute the difference and find the sum of coefficients
    difference = V1_special - V1_SPE
    
    # Simplify and get the numerator and denominator
    difference = sympy.simplify(difference)
    p_d, q_d = sympy.fraction(difference)

    # Ensure p(d) and q(d) have integer coefficients by clearing denominators if any
    p_d, q_d = sympy.poly(p_d, d).LC() * p_d, sympy.poly(q_d, d).LC() * q_d
    p_d, q_d = sympy.fraction(sympy.simplify(p_d / q_d))

    # To get a standard form, let's use the original derivation before simplification
    # p(d) = (4 - 3d)(1 + d) - (4 - 2d)
    # q(d) = (4 - 2d)(1 + d)
    p_d_unsimplified = (4 - 3*d)*(1 + d) - 1*(4 - 2*d)
    q_d_unsimplified = (4 - 2*d)*(1 + d)
    
    p_d_expanded = sympy.expand(p_d_unsimplified)
    q_d_expanded = sympy.expand(q_d_unsimplified)
    
    # Extract coefficients
    p_coeffs = sympy.Poly(p_d_expanded, d).all_coeffs()
    q_coeffs = sympy.Poly(q_d_expanded, d).all_coeffs()
    
    # Sum of absolute values of coefficients
    sum_abs_coeffs_p = sum(abs(c) for c in p_coeffs)
    sum_abs_coeffs_q = sum(abs(c) for c in q_coeffs)
    total_sum = sum_abs_coeffs_p + sum_abs_coeffs_q
    
    print("Player 1's payoff in SPE: 1/(1+d)")
    print(f"Player 2's simulated payoff from rejection: {rejection_threshold}")
    print(f"Player 1's optimal payoff against specific P2: {V1_special}")
    print("\n--- Difference Calculation ---")
    print(f"Difference D(d) = ({p_d_expanded}) / ({q_d_expanded})")
    
    print(f"\nNumerator p(d) = {p_d_expanded}")
    p_coeffs_int = [int(c) for c in p_coeffs]
    print(f"Numerator coefficients: {p_coeffs_int}")
    print(f"Sum of absolute values of numerator coefficients: {int(sum_abs_coeffs_p)}")

    print(f"\nDenominator q(d) = {q_d_expanded}")
    q_coeffs_int = [int(c) for c in q_coeffs]
    print(f"Denominator coefficients: {q_coeffs_int}")
    print(f"Sum of absolute values of denominator coefficients: {int(sum_abs_coeffs_q)}")
    
    print(f"\nTotal sum of absolute values of coefficients: {int(total_sum)}")

solve_bargaining_problem()