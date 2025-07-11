import sympy

def solve_bargaining_problem():
    """
    Solves the described bargaining problem by performing symbolic calculations.
    """
    # Define the discount factor 'd' as a symbolic variable
    d = sympy.Symbol('d')

    # 1. Player 1's payoff in the standard SPE
    V_SPE = 1 / (1 + d)
    print(f"Step 1: Player 1's payoff in the standard SPE is V_SPE = {V_SPE}\n")

    # 2. Analyze Player 2's reservation value based on the R-simulation
    # Let F1 be P2's expected payoff share when P1 proposes in the simulation.
    # Let F2 be P2's expected payoff share when P2 proposes in the simulation.
    # In the simulation, offers are U([0,1]) (so E[offer]=1/2) and A/R is 50/50.
    # F1 = (1/2)*E[1-z] + (1/2)*d*F2 = (1/2)*(1/2) + (d/2)*F2 = 1/4 + d/2*F2
    # F2 = (1/2)*E[1-y] + (1/2)*d*F1 = (1/2)*(1/2) + (d/2)*F1 = 1/4 + d/2*F1
    F1, F2 = sympy.symbols('F1 F2')
    eq1 = sympy.Eq(F1, sympy.Rational(1, 4) + (d/2)*F2)
    eq2 = sympy.Eq(F2, sympy.Rational(1, 4) + (d/2)*F1)

    # Solve the system for F2
    solution = sympy.solve([eq1, eq2], [F1, F2])
    F2_val = solution[F2]
    
    # P2's reservation payoff comes from rejecting P1's offer in Round 1.
    # The pie becomes 'd' and it's P2's turn to propose in the simulation.
    # So the reservation payoff is d * F2_val
    P2_reservation_payoff = d * F2_val
    P2_reservation_payoff_simplified = sympy.simplify(P2_reservation_payoff)
    print(f"Step 2: Player 2's simulated payoff share (F2) is {F2_val}")
    print(f"Step 3: Player 2's reservation payoff in the actual game is {P2_reservation_payoff_simplified}\n")

    # 4. Player 1's optimal payoff
    # P1 offers a split (x, 1-x) where 1-x equals P2's reservation value.
    # P1's payoff is x = 1 - P2_reservation_payoff
    V_optimal = 1 - P2_reservation_payoff
    V_optimal_simplified = sympy.simplify(V_optimal)
    print(f"Step 4: Player 1's optimal payoff is V_optimal = 1 - {P2_reservation_payoff_simplified} = {V_optimal_simplified}\n")

    # 5. Calculate the difference
    difference = V_optimal - V_SPE
    # Use cancel to get the simplified rational function form p(d)/q(d)
    difference_simplified = sympy.cancel(difference)
    print(f"Step 5: The difference V_optimal - V_SPE is {difference_simplified}\n")

    # 6. Find the sum of absolute values of coefficients
    p, q = sympy.fraction(difference_simplified)

    # For a canonical representation, make the leading coefficient of the denominator positive
    if q.as_poly(d).LC() < 0:
        p = -p
        q = -q
        
    p_poly = sympy.Poly(p, d)
    q_poly = sympy.Poly(q, d)

    p_coeffs = [int(c) for c in p_poly.all_coeffs()]
    q_coeffs = [int(c) for c in q_poly.all_coeffs()]

    print("Step 6: Expressing the difference as p(d)/q(d):")
    # Using as_expr() to display the polynomial in a standard way
    print(f"The numerator p(d) is: {p_poly.as_expr()}")
    print(f"The denominator q(d) is: {q_poly.as_expr()}")
    
    print("\nThe integer coefficients for the final equation are:")
    print(f"Numerator coefficients: {p_coeffs}")
    print(f"Denominator coefficients: {q_coeffs}")

    p_sum_abs = sum(abs(c) for c in p_coeffs)
    q_sum_abs = sum(abs(c) for c in q_coeffs)
    total_sum_abs = p_sum_abs + q_sum_abs

    print(f"\nSum of absolute values of coefficients:")
    print(f"Numerator: {p_sum_abs}")
    print(f"Denominator: {q_sum_abs}")
    print(f"Total Sum: {total_sum_abs}")

solve_bargaining_problem()
<<<M>>>