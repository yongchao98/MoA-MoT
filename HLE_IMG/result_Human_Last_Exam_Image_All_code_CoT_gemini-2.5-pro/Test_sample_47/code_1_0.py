import numpy as np

def solve_circuit():
    """
    This function calculates the current I in the given electrical circuit.
    """
    # --- Define component values from the diagram ---
    V_source = 41.0  # Voltage in Volts
    R76 = 76.0
    R8 = 8.0
    R14 = 14.0
    R11 = 11.0
    R29 = 29.0

    # --- Analysis Steps ---

    # Step 1: Calculate the equivalent resistance of the parallel resistors R76 and R8.
    # These resistors are connected between the 41V source and a node we'll call Vx.
    R_p = (R76 * R8) / (R76 + R8)

    # Step 2: Set up the system of linear equations for nodal analysis.
    # We have two unknown node voltages: Vx and Vy.
    # Let Vx be the voltage at the node connecting R_p and R14.
    # Let Vy be the voltage at the node connecting R14, R11, and R29.
    # The system is in the form A * V = B, where V = [Vx, Vy]'

    # KCL equation at Node Vx: (Vx - V_source)/R_p + (Vx - Vy)/R14 = 0
    # KCL equation at Node Vy: (Vy - Vx)/R14 + (Vy - V_source)/R11 + Vy/R29 = 0
    # Rearranging into the form A*V = B:
    # Eq1: Vx*(1/R_p + 1/R14) - Vy*(1/R14) = V_source/R_p
    # Eq2: -Vx*(1/R14) + Vy*(1/R14 + 1/R11 + 1/R29) = V_source/R11

    # Coefficient matrix A
    A = np.array([
        [1/R_p + 1/R14, -1/R14],
        [-1/R14, 1/R14 + 1/R11 + 1/R29]
    ])

    # Constant vector B
    B = np.array([
        V_source / R_p,
        V_source / R11
    ])

    # Step 3: Solve for the node voltages Vx and Vy using numpy.
    try:
        voltages = np.linalg.solve(A, B)
        Vx = voltages[0]
        Vy = voltages[1]

        # Step 4: Calculate the current I.
        # I is the current flowing from Node X to Node Y through R14.
        I = (Vx - Vy) / R14

        # --- Output the results ---
        print("Step 1: Calculate equivalent resistance of parallel resistors (76 Ohm and 8 Ohm).")
        print(f"R_parallel = ({R76} * {R8}) / ({R76} + {R8}) = {R_p:.4f} Ohms.")
        
        print("\nStep 2: Solve for unknown node voltages Vx and Vy using Nodal Analysis.")
        print(f"The calculated node voltages are:")
        print(f"Vx = {Vx:.4f} V")
        print(f"Vy = {Vy:.4f} V")

        print("\nStep 3: Calculate the current I using Ohm's Law.")
        print("I is the current flowing through the 14 Ohm resistor.")
        print(f"I = (Vx - Vy) / R14")
        print(f"I = ({Vx:.4f} V - {Vy:.4f} V) / {R14} Ohm")
        print(f"I = {I:.4f} A")

    except np.linalg.LinAlgError:
        print("The system of equations is singular, and the circuit cannot be solved as configured.")

# Run the solver
solve_circuit()