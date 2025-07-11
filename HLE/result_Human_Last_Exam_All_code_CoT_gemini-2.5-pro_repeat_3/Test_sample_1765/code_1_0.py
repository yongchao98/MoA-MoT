import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal QSH device
    using the Landauer-Büttiker formalism.
    """
    print("This script calculates the two-terminal conductance of a four-terminal device")
    print("exhibiting the Quantum Spin Hall effect with helical edge states.")
    print("-" * 60)
    
    # 1. & 2. Define symbols and Landauer-Büttiker equations
    # V_i are the voltages at the terminals.
    # I is the injection current.
    # G0 represents the quantum of conductance, e^2/h.
    V1, V2, V3, V4, I, G0 = sympy.symbols('V1, V2, V3, V4, I, G0')

    print("The Landauer-Büttiker equations for the system are:")
    # For each terminal 'i', the current is I_i = G0 * (N_i*V_i - sum(T_ij*V_j)),
    # where N_i is the number of channels leaving terminal i (here N_i=2 for all i)
    # and T_ij is the transmission from terminal j to i.
    # Based on the helical edge states (spin-up CW, spin-down CCW):
    # - Channels to 1: from 2 (spin-down) and 4 (spin-up)
    # - Channels to 2: from 1 (spin-up) and 3 (spin-down)
    # - Channels to 3: from 2 (spin-up) and 4 (spin-down)
    # - Channels to 4: from 1 (spin-down) and 3 (spin-up)
    
    eq1 = sympy.Eq(G0 * (2*V1 - V2 - V4), I)    # I_1 = I
    eq2 = sympy.Eq(G0 * (2*V2 - V1 - V3), -I)   # I_2 = -I
    eq3 = sympy.Eq(G0 * (2*V3 - V2 - V4), 0)    # I_3 = 0 (floated)
    eq4 = sympy.Eq(G0 * (2*V4 - V1 - V3), 0)    # I_4 = 0 (floated)
    
    print(f"I_1 = G0*(2*V1 - V2 - V4) = I")
    print(f"I_2 = G0*(2*V2 - V1 - V3) = -I")
    print(f"I_3 = G0*(2*V3 - V2 - V4) = 0")
    print(f"I_4 = G0*(2*V4 - V1 - V3) = 0")
    print("-" * 60)

    # 3. Solve the system of equations.
    # We want to find the relationship between I and (V1 - V2).
    # We solve for the potentials V1, V2, V3, and the current I in terms of one potential, say V4.
    # The system is linear, so sympy can solve it.
    solution = sympy.solve([eq1, eq2, eq3, eq4], (V1, V2, V3, I))

    # The solution will be a dictionary, e.g., {V1: expr_for_V1, V2: expr_for_V2, ...}
    # All expressions will be in terms of V4 and G0.
    v1_expr = solution[V1]
    v2_expr = solution[V2]
    i_expr = solution[I]
    
    # 4. Calculate the conductance G_12 = I / (V1 - V2)
    voltage_difference = v1_expr - v2_expr
    conductance = i_expr / voltage_difference
    
    # The result will be in terms of G0. Sympy simplifies the expression.
    final_conductance = sympy.simplify(conductance)
    
    # 5. Print the final result in the requested format.
    print("By solving this system, we find the conductance G_12 = I / (V1 - V2).")
    
    num, den = sympy.fraction(final_conductance.as_expr()/G0)
    
    print("\nThe final equation for the conductance is:")
    print("G_12 = ", end="")
    print(num, end="")
    print(" / ", end="")
    print(den, end="")
    print(" * e^2/h")

calculate_qsh_conductance()