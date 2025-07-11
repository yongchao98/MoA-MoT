import numpy as np
from fractions import Fraction

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device with two terminals floated.
    """
    
    # 1. Define the experimental conditions.
    # We apply a voltage V across terminals 1 and 2. For simplicity,
    # let's set V1 = 1 and V2 = 0. The conductance is independent of V.
    V1 = 1.0
    V2 = 0.0
    
    # 2. Set up the system of linear equations for the floating terminals.
    # The condition that terminals 3 and 4 are floating (I3=0, I4=0) gives:
    # -2*V3 + 1*V4 = -V2
    #  1*V3 - 2*V4 = -V1
    # We can write this in matrix form A*x = b, where x = [V3, V4].
    
    A = np.array([
        [-2,  1],
        [ 1, -2]
    ])
    
    b = np.array([
        -V2,
        -V1
    ])
    
    # 3. Solve for the unknown voltages V3 and V4.
    try:
        floating_voltages = np.linalg.solve(A, b)
        V3 = floating_voltages[0]
        V4 = floating_voltages[1]
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")
        return

    # 4. Calculate the current I flowing from terminal 1 to 2.
    # The current is I = I2 = G0 * (V1 + V3 - 2*V2), where G0 = e^2/h.
    # We calculate the numerical coefficient for G0.
    current_coefficient = V1 + V3 - (2 * V2)
    
    # 5. Calculate the two-terminal conductance G12 = I / (V1 - V2).
    # G12 = (current_coefficient * G0) / (V1 - V2)
    conductance_coefficient = current_coefficient / (V1 - V2)
    
    # 6. Print the final result in a clear format.
    # We use the Fraction class to get a nice rational number.
    conductance_fraction = Fraction(conductance_coefficient).limit_denominator()
    numerator = conductance_fraction.numerator
    denominator = conductance_fraction.denominator
    
    print("Based on the Landauer-Büttiker formalism for a four-terminal device with helical edge states:")
    print("The two-terminal conductance G_12 between terminals 1 and 2 with terminals 3 and 4 floated is:")
    print(f"G_12 = ({numerator}/{denominator}) * e^2/h")
    print("\nIn this final equation:")
    print(f"The first number is the numerator: {numerator}")
    print(f"The second number is the denominator: {denominator}")

    # Return the final numerical answer for the platform.
    global final_answer
    final_answer = conductance_coefficient


# Execute the calculation and store the result
final_answer = None
calculate_qsh_conductance()

# The final answer format required by the platform.
# Using print() to make sure it's part of the script's output.
print(f"\n<<<final_answer>>>")
print(f"<<<{final_answer}>>>")