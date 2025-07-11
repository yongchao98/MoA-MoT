import numpy as np

def solve_circuit():
    """
    Solves the given electrical circuit to find the current I.
    """
    # Step 1: Define circuit parameters
    V_source = 41  # V
    R76 = 76       # Ohms
    R8 = 8         # Ohms
    R14 = 14       # Ohms
    R11 = 11       # Ohms
    R29 = 29       # Ohms

    # Step 2: Calculate the equivalent resistance of the parallel resistors on the left
    R_L_equiv = (R76 * R8) / (R76 + R8)

    # Step 3: Set up the system of linear equations for nodal analysis (A*x = b)
    # The unknowns are x = [V_L, V_R]
    # KCL at Node L: (V_source - V_L)/R_L_equiv + (V_R - V_L)/R14 = 0
    # KCL at Node R: (V_source - V_R)/R11 + (V_L - V_R)/R14 + (0 - V_R)/R29 = 0
    # Rearranging into matrix form:
    # Eq1: V_L * (1/R_L_equiv + 1/R14) - V_R * (1/R14) = V_source / R_L_equiv
    # Eq2: -V_L * (1/R14) + V_R * (1/R11 + 1/R14 + 1/R29) = V_source / R11
    
    A = np.array([
        [1/R_L_equiv + 1/R14, -1/R14],
        [-1/R14, 1/R11 + 1/R14 + 1/R29]
    ])
    
    b = np.array([
        V_source / R_L_equiv,
        V_source / R11
    ])

    # Step 4: Solve for the node voltages V_L and V_R
    try:
        voltages = np.linalg.solve(A, b)
        V_L = voltages[0]
        V_R = voltages[1]
    except np.linalg.LinAlgError:
        print("The system of equations is singular. Cannot solve.")
        return

    # Step 5: Calculate the current I
    I = (V_L - V_R) / R14

    # Print the final equation with all the numbers
    print("The value of the current I is found using the node voltages V_L and V_R and Ohm's law.")
    print(f"V_L = {V_L:.2f} V")
    print(f"V_R = {V_R:.2f} V")
    print("I = (V_L - V_R) / R14")
    print(f"I = ({V_L:.2f} - {V_R:.2f}) / {R14} = {I:.3f} A")

solve_circuit()