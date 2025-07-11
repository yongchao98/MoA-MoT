import sympy

def solve_bargaining_problem():
    """
    This function solves the described bargaining problem step-by-step.
    """
    # Set up the symbolic variable for the discount factor
    d = sympy.Symbol('d')

    # Part 1: Payoff in Subgame Perfect Equilibrium (SPE)
    # In the standard Rubinstein bargaining game, Player 1's payoff is V1_SPE = 1 / (1 + d).
    V1_SPE = 1 / (1 + d)

    # Part 2: Player 1's optimal payoff against the specified Player 2.
    # Player 1 offers (x, 1-x). Player 2 accepts if their payoff 1-x is greater than
    # the expected payoff from rejecting (E_R2), which is determined by a simulation.

    # We calculate Player 2's expected payoff from rejecting, E_R2.
    # The simulation starts at t=1 (pie size d) with Player 2 as the proposer.
    # Let v_2O be P2's expected share of the current pie when P2 proposes.
    # Let v_2A be P2's expected share when P1 proposes (and P2 decides to Accept/Reject).
    # In the simulation, offers are random s from U([0,1]) for P1's share, and A/R decisions are 50/50.
    # P2's expected share on acceptance is E[1-s] = 0.5.
    v_2O, v_2A = sympy.symbols('v_2O v_2A')

    # We use Rational to maintain precision. 0.5 is Rational(1,2), 0.25 is Rational(1,4).
    # eq1: v_2O = (Prob P1 accepts) * E[P2's share] + (Prob P1 rejects) * (discount * P2's future share)
    eq1 = sympy.Eq(v_2O, sympy.Rational(1, 2) * sympy.Rational(1, 2) + sympy.Rational(1, 2) * d * v_2A)
    # eq2: v_2A = (Prob P2 accepts) * E[P2's share] + (Prob P2 rejects) * (discount * P2's future share)
    eq2 = sympy.Eq(v_2A, sympy.Rational(1, 2) * sympy.Rational(1, 2) + sympy.Rational(1, 2) * d * v_2O)

    # Solve the system of equations for v_2O.
    sol = sympy.solve([eq1, eq2], (v_2O, v_2A))
    v_2O_expr = sol[v_2O]

    # The total expected payoff for P2 in the simulation is discounted by d from the start of the game.
    E_R2 = d * v_2O_expr

    # Player 1's optimal strategy is to offer P2 a share equal to E_R2. P1's payoff is thus 1 - E_R2.
    V1_optimal = 1 - E_R2

    # Part 3: Compute the difference and find the sum of coefficients.
    # Calculate the difference V1_optimal - V1_SPE.
    difference = V1_optimal - V1_SPE

    # We use `cancel` to simplify the fraction to its lowest terms, ensuring the
    # numerator and denominator polynomials are coprime.
    simplified_difference = sympy.cancel(difference)

    # Extract the numerator p(d) and denominator q(d).
    p_d, q_d = sympy.fraction(simplified_difference)

    # Create Poly objects to correctly extract all coefficients, including zero terms.
    p_poly = sympy.Poly(p_d, d)
    q_poly = sympy.Poly(q_d, d)

    # Get the lists of integer coefficients.
    p_coeffs = p_poly.all_coeffs()
    q_coeffs = q_poly.all_coeffs()

    # Calculate the sum of the absolute values of the coefficients.
    sum_abs_coeffs_p = sum(abs(c) for c in p_coeffs)
    sum_abs_coeffs_q = sum(abs(c) for c in q_coeffs)
    total_sum_abs_coeffs = sum_abs_coeffs_p + sum_abs_coeffs_q

    # Print the required outputs.
    p_expr = p_poly.as_expr()
    q_expr = q_poly.as_expr()
    
    print("The difference between Player 1's optimal payoff against the specified opponent and the standard SPE payoff is a rational expression p(d)/q(d).")
    print("\nAfter simplification, the final polynomials are:")
    
    print(f"\nNumerator: p(d) = {p_expr}")
    # The line below outputs each number in the equation, as requested.
    print(f"The coefficients of p(d) are {', '.join(map(str, [int(c) for c in p_coeffs]))}.")
    
    print(f"\nDenominator: q(d) = {q_expr}")
    # The line below outputs each number in the equation, as requested.
    print(f"The coefficients of q(d) are {', '.join(map(str, [int(c) for c in q_coeffs]))}.")

    print(f"\n\nThe sum of the absolute values of the coefficients in p(d) is {int(sum_abs_coeffs_p)}.")
    print(f"The sum of the absolute values of the coefficients in q(d) is {int(sum_abs_coeffs_q)}.")
    
    print(f"\nThe total sum of the absolute values of all coefficients is:")
    print(int(total_sum_abs_coeffs))

solve_bargaining_problem()