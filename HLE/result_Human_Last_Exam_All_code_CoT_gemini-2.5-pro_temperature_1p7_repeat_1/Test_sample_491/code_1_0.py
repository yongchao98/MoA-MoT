import sympy

def solve_bargaining_problem():
    """
    Solves the bargaining problem as described.
    1.  Calculates the SPE payoff for Player 1.
    2.  Calculates the optimal payoff for Player 1 against the specified opponent.
    3.  Finds the difference and expresses it as p(d)/q(d).
    4.  Sums the absolute values of the integer coefficients of p(d) and q(d).
    """
    # Define the symbolic variable for the discount factor
    d = sympy.Symbol('d')

    # Part 1: Payoff in the standard Subgame Perfect Equilibrium (SPE)
    # In the standard Rubinstein-Stahl model, Player 1's payoff is V1 = 1 / (1 + d).
    V1_SPE = 1 / (1 + d)

    # Part 2: Player 1's optimal payoff against the specific opponent.
    # Player 2's decision to accept/reject depends on a simulated payoff from rejecting.
    # The simulated payoff, U_R, is based on a game with random actions.
    # Let E22 be the expected payoff for Player 2, in a random simulation,
    # when it is Player 2's turn to propose with a pie of size 1.
    # This value is derived by solving the system of equations for the random game:
    # E2_2 = 1/4 + (d/2)*E2_1 and E2_1 = 1/4 + (d/2)*E2_2
    # The solution is E2_2 = 1 / (4 - 2*d).
    E22 = 1 / (4 - 2 * d)

    # The simulated payoff for Player 2 upon rejecting Player 1's initial offer is U_R = d * E22,
    # as the pie has been discounted to size d.
    U_R = d * E22

    # Player 2 will accept an offer (x, 1-x) if their real payoff 1-x is greater than
    # their simulated rejection payoff U_R.
    # Player 1's optimal strategy is to make the highest possible offer x that P2 will accept,
    # which is when 1-x = U_R. So, V1_optimal = x = 1 - U_R.
    V1_optimal = 1 - U_R

    # Part 3: Compute the difference
    # Difference = V1_optimal - V1_SPE
    Difference = V1_optimal - V1_SPE

    # Part 4: Express as a quotient of polynomials p(d)/q(d)
    # sympy.fraction() finds the numerator and denominator over a common base.
    p_d, q_d = sympy.fraction(Difference)

    # Expand the polynomials to work with their standard form.
    p_d_expanded = sympy.expand(p_d)
    q_d_expanded = sympy.expand(q_d)

    # Part 5: Sum the absolute values of the coefficients
    p_poly = sympy.Poly(p_d_expanded, d)
    q_poly = sympy.Poly(q_d_expanded, d)

    p_coeffs = p_poly.all_coeffs()
    q_coeffs = q_poly.all_coeffs()

    sum_abs_coeffs = sum(abs(c) for c in p_coeffs) + sum(abs(c) for c in q_coeffs)

    # Print the breakdown of the final equation and the result as requested
    p_abs_coeffs_str = [f"|{int(c)}|" for c in p_coeffs]
    q_abs_coeffs_str = [f"|{int(c)}|" for c in q_coeffs]

    print("The difference between the two payoffs is p(d)/q(d).")
    print(f"Numerator: p(d) = {p_d_expanded}")
    print(f"Denominator: q(d) = {q_d_expanded}")
    print("\nThe sum of the absolute values of the integer coefficients is calculated as follows:")
    final_equation_str = " + ".join(p_abs_coeffs_str) + " + " + " + ".join(q_abs_coeffs_str)
    print(f"{final_equation_str} = {int(sum_abs_coeffs)}")

if __name__ == '__main__':
    solve_bargaining_problem()