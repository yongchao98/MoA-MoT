import sympy

def solve_bargaining_problem():
    """
    Solves the specified bargaining problem by calculating two different payoffs for Player 1,
    finding their difference, and summing the absolute values of the coefficients of the
    resulting rational polynomial expression.
    """
    # Define the discount factor 'd' as a symbolic variable
    d = sympy.Symbol('d')

    # Part 1: Calculate Player 1's payoff in the unique subgame perfect equilibrium (SPE)
    V1, V2 = sympy.symbols('V1 V2')
    spe_eq1 = sympy.Eq(V1, 1 - d * V2)
    spe_eq2 = sympy.Eq(V2, 1 - d * V1)
    spe_solution = sympy.solve([spe_eq1, spe_eq2], (V1, V2))
    p1_spe_payoff = spe_solution[V1]

    # Part 2: Calculate Player 1's optimal payoff against the specific Player 2
    
    # First, determine the payoff Player 2 expects from the 'Reject' simulation.
    # E2p: P2's expected payoff in sim when proposing
    # E2r: P2's expected payoff in sim when responding
    E2p, E2r = sympy.symbols('E2p E2r')
    
    # In a sim, offer is from U([0,1]), so expected share is 0.5. Other's share is 1-0.5=0.5.
    # A/R decisions are 50/50.
    sim_eq1 = sympy.Eq(E2p, sympy.Rational(1, 2) * sympy.Rational(1, 2) + sympy.Rational(1, 2) * d * E2r)
    sim_eq2 = sympy.Eq(E2r, sympy.Rational(1, 2) * sympy.Rational(1, 2) + sympy.Rational(1, 2) * d * E2p)

    sim_solution = sympy.solve([sim_eq1, sim_eq2], (E2p, E2r))
    e2p_value = sim_solution[E2p]

    # The payoff for P2 in the 'Reject' simulation is the discounted value of proposing.
    p2_rejection_sim_payoff = d * e2p_value
    
    # Player 1 offers Player 2 exactly this amount to ensure acceptance.
    # Player 1's payoff is 1 minus what is offered to Player 2.
    p1_optimal_payoff = 1 - p2_rejection_sim_payoff
    
    # Part 3: Compute the difference and analyze the coefficients.
    
    # Calculate the difference between the optimal payoff and the SPE payoff
    difference = p1_optimal_payoff - p1_spe_payoff
    
    # Simplify the expression to a canonical p(d)/q(d) form
    simplified_difference = sympy.cancel(difference)
    
    # Extract the numerator and denominator
    p_d, q_d = simplified_difference.as_numer_denom()

    # To ensure a standard form, make the leading coefficient of the denominator positive
    if q_d.as_poly(d).LC() < 0:
        p_d = -p_d
        q_d = -q_d

    # Get the coefficients from the polynomials
    p_coeffs = p_d.as_poly(d).all_coeffs()
    q_coeffs = q_d.as_poly(d).all_coeffs()

    # Sum the absolute values of the coefficients
    sum_abs_coeffs = sum(abs(c) for c in p_coeffs) + sum(abs(c) for c in q_coeffs)
    
    # Output the steps as per the instructions
    print(f"The final expression for the difference is p(d)/q(d), where:")
    print(f"p(d) = {p_d}")
    print(f"q(d) = {q_d}")
    print(f"\nBreaking this down into coefficients for the final sum:")
    print(f"Coefficients of the numerator p(d) are: {p_coeffs}")
    print(f"Coefficients of the denominator q(d) are: {q_coeffs}")
    
    # Print each number in the final equation's coefficient sums
    abs_p_coeffs = [abs(c) for c in p_coeffs]
    abs_q_coeffs = [abs(c) for c in q_coeffs]
    p_sum_str = " + ".join(map(str, abs_p_coeffs))
    q_sum_str = " + ".join(map(str, abs_q_coeffs))
    
    print(f"The sum is ({p_sum_str}) + ({q_sum_str})")
    
    print(f"\nThe total sum of the absolute values of the coefficients is: {sum_abs_coeffs}")

solve_bargaining_problem()