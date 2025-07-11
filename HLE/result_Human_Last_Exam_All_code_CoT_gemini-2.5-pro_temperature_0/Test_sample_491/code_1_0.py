import sympy

def solve_bargaining_problem():
    """
    Solves the described bargaining problem by:
    1. Defining the SPE and modified game payoffs for Player 1.
    2. Calculating the difference between them.
    3. Expressing the difference as a quotient of polynomials p(d)/q(d).
    4. Calculating the sum of the absolute values of the coefficients of p(d) and q(d).
    """
    # Define the discount factor 'd' as a symbolic variable
    d = sympy.Symbol('d')

    # Part (2): Payoff in the unique subgame perfect equilibrium (SPE)
    V_spe = 1 / (1 + d)

    # Part (1): Payoff for Player 1 against the specific opponent
    # P2's decision is based on a simulation where expected payoff for the proposer is E_p = 1/(2*(2-d)).
    # The R-simulation starts with P2 proposing on a pie of size d.
    # So, P2's expected payoff from rejecting is d * E_p.
    E_R_sim = d / (2 * (2 - d))
    
    # P1 offers P2 just enough to accept, so 1 - V_mod = E_R_sim
    V_mod = 1 - E_R_sim
    
    # Simplify V_mod
    V_mod = sympy.simplify(V_mod)
    
    # Calculate the difference D = V_mod - V_spe
    difference = V_mod - V_spe

    # The problem asks to express the difference as p(d)/q(d).
    # sympy.cancel() simplifies the fraction and returns it in a standard form.
    simplified_difference = sympy.cancel(difference)

    # Extract the numerator p(d) and denominator q(d)
    p_d, q_d = sympy.fraction(simplified_difference)

    # Convert them to Poly objects to easily access coefficients
    p_poly = sympy.Poly(p_d, d)
    q_poly = sympy.Poly(q_d, d)

    # Get all coefficients
    p_coeffs = p_poly.all_coeffs()
    q_coeffs = q_poly.all_coeffs()

    # Calculate the sum of the absolute values of the coefficients
    p_abs_sum = sum(abs(c) for c in p_coeffs)
    q_abs_sum = sum(abs(c) for c in q_coeffs)
    total_sum = p_abs_sum + q_abs_sum

    # Print the results step-by-step as requested
    print(f"Player 1's payoff in SPE: V_spe = {V_spe}")
    print(f"Player 1's payoff against the specific opponent: V_mod = {V_mod}")
    print(f"The difference is D = V_mod - V_spe = {simplified_difference}")
    print("-" * 20)
    print(f"Numerator p(d) = {p_d}")
    print(f"Coefficients of p(d): {p_coeffs}")
    print(f"Sum of absolute values of coefficients in p(d): {p_abs_sum}")
    print("-" * 20)
    print(f"Denominator q(d) = {q_d}")
    print(f"Coefficients of q(d): {q_coeffs}")
    print(f"Sum of absolute values of coefficients in q(d): {q_abs_sum}")
    print("-" * 20)
    print(f"Total sum of absolute values of coefficients = {p_abs_sum} + {q_abs_sum} = {total_sum}")

solve_bargaining_problem()