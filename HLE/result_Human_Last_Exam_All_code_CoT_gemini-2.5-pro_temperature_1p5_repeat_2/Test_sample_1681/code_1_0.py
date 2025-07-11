import sympy

def solve_quantum_prisoners_dilemma():
    """
    Calculates the optimal equilibrium point for the quantum prisoner's dilemma.

    This function follows the EWL protocol to find the condition for a new
    Nash Equilibrium (Q, Q) to exist and calculates the payoff at the optimal
    point of that equilibrium.
    """
    # Step 1: Define payoffs from the matrix
    # R: (Cooperate, Cooperate), S: (Cooperate, Defect)
    # T: (Defect, Cooperate), P: (Defect, Defect)
    R, S, T, P = 5, 0, 7, 1

    # Step 2: Define symbolic variables for entanglement
    # sg2 represents sin^2(gamma/2)
    # cg2 represents cos^2(gamma/2)
    sg2, cg2 = sympy.symbols('sg2 cg2')

    # Step 3: Formulate payoff equations based on the EWL protocol
    # Payoff for player A if both play the quantum strategy 'Q'
    payoff_QQ = R * cg2 + P * sg2
    
    # Payoff for player A if she defects ('D') while player B plays 'Q'
    payoff_DQ = T * sg2 + S * cg2

    # Step 4: Establish the Nash Equilibrium condition
    # The payoff for playing Q must be at least as good as defecting
    # A(Q, Q) >= A(D, Q)
    equilibrium_condition = sympy.Ge(payoff_QQ, payoff_DQ)

    # Step 5: Solve for the optimal entanglement
    # Substitute cg2 = 1 - sg2 to solve in terms of sg2
    condition_in_sg2 = equilibrium_condition.subs(cg2, 1 - sg2)
    
    # The inequality becomes: 5*(1-sg2) + 1*sg2 >= 7*sg2
    # 5 - 5*sg2 + sg2 >= 7*sg2
    # 5 - 4*sg2 >= 7*sg2
    # 5 >= 11*sg2  => sg2 <= 5/11
    
    # We choose the boundary condition to maximize the equilibrium payoff
    optimal_sg2 = sympy.Rational(5, 11)
    optimal_cg2 = 1 - optimal_sg2

    # Step 6: Calculate the final equilibrium payoff
    # Substitute the optimal entanglement value back into the payoff equation for (Q, Q)
    equilibrium_payoff_expr = R * cg2 + P * sg2
    equilibrium_payoff = equilibrium_payoff_expr.subs({cg2: optimal_cg2, sg2: optimal_sg2})

    # Display the final calculation steps
    # The payoff formula simplifies to 4*cg2 + 1
    # Payoff = 4 * (6/11) + 1 = 24/11 + 11/11 = 35/11
    p1 = 4
    p2 = optimal_cg2.p
    p3 = optimal_cg2.q
    p4 = 1
    
    print("The equilibrium is found by choosing an optimal level of entanglement.")
    print("The condition for the quantum strategy (Q,Q) to be a Nash Equilibrium is sin^2(gamma/2) <= 5/11.")
    print("Choosing the maximum allowed entanglement, sin^2(gamma/2) = 5/11, we get cos^2(gamma/2) = 6/11.")
    print("\nThe payoff at this equilibrium point is calculated using the formula: R*cos^2(gamma/2) + P*sin^2(gamma/2)")
    print("Which simplifies to: 4*cos^2(gamma/2) + 1")
    print("\nFinal Calculation:")
    print(f"{p1} * ({p2}/{p3}) + {p4} = {equilibrium_payoff}")

solve_quantum_prisoners_dilemma()
<<<35/11>>>