import numpy as np

def solve_circuit():
    """
    Analyzes the given electrical circuit to find the current I.
    """
    # Given values from the circuit diagram
    V_source = 41.0  # Voltage in Volts
    R76 = 76.0       # Resistance in Ohms
    R8 = 8.0         # Resistance in Ohms
    R14 = 14.0       # Resistance in Ohms
    R11 = 11.0       # Resistance in Ohms
    R29 = 29.0       # Resistance in Ohms

    # Step 1: Calculate the equivalent resistance of the parallel resistors R76 and R8.
    # This combination is connected between the 41V source and Node 1.
    R_p1 = (R76 * R8) / (R76 + R8)

    # Step 2: Set up the matrix equations for nodal analysis (A * V = b),
    # where V = [V1, V2] are the voltages at the two unknown nodes.
    # KCL at Node 1: (V1 - V_source)/R_p1 + (V1 - V2)/R14 = 0
    # => V1*(1/R_p1 + 1/R14) - V2*(1/R14) = V_source/R_p1
    #
    # KCL at Node 2: (V2 - V1)/R14 + (V2 - V_source)/R11 + V2/R29 = 0
    # => -V1*(1/R14) + V2*(1/R14 + 1/R11 + 1/R29) = V_source/R11

    # Matrix A contains the coefficients of V1 and V2
    A = np.array([
        [1/R_p1 + 1/R14, -1/R14],
        [-1/R14, 1/R14 + 1/R11 + 1/R29]
    ])

    # Vector b contains the constant terms of the equations
    b = np.array([
        V_source / R_p1,
        V_source / R11
    ])

    try:
        # Step 3: Solve the system of linear equations for the node voltages V1 and V2
        node_voltages = np.linalg.solve(A, b)
        V1 = node_voltages[0]
        V2 = node_voltages[1]

        # Step 4: Calculate the current I flowing from Node 1 to Node 2 through R14
        I = (V1 - V2) / R14

        # Print the results and the final equation
        print("Step 1: Find the equivalent resistance of the parallel 76 Ohm and 8 Ohm resistors.")
        print(f"R_p1 = ({R76} * {R8}) / ({R76} + {R8}) = {R_p1:.4f} Ohms")
        
        print("\nStep 2: Solve for the node voltages V1 and V2 using nodal analysis.")
        print(f"The calculated node voltages are:")
        print(f"V1 = {V1:.4f} V")
        print(f"V2 = {V2:.4f} V")

        print("\nStep 3: Calculate the current I, which flows from Node 1 to Node 2 through the 14 Ohm resistor.")
        print(f"I = (V1 - V2) / R14")
        print(f"I = ({V1:.4f} - {V2:.4f}) / {R14}")
        print(f"I = {I:.4f} A")

    except np.linalg.LinAlgError:
        print("The system of equations is singular and cannot be solved.")

solve_circuit()