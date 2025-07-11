import sys

def solve_circuit():
    """
    Calculates the total current flowing from the voltage source in the given circuit.
    """
    # Define the component values from the circuit diagram
    V2 = 1.0  # Voltage in Volts
    R1 = 3.0  # Resistance in Ohms
    R2 = 7.0  # Resistance in Ohms
    R3 = 9.0  # Resistance in Ohms
    R7 = 100.0 # Resistance in Ohms
    R8 = 100.0 # Resistance in Ohms

    # The circuit has two main parallel branches. The total current is the sum of the currents in these branches.
    # The final equation for the total current (I) is:
    # I = I_branch1 + I_branch2
    # I = V2/R3 + V2 / ( (1/(1/R1 + 1/R2 + 1/R8)) + R7 )

    # Calculate the current in the first branch (containing only R3)
    I_branch1 = V2 / R3

    # Calculate the equivalent resistance of the parallel combination of R1, R2, and R8
    R_p128 = 1 / (1/R1 + 1/R2 + 1/R8)
    
    # Calculate the total resistance of the second branch (R_p128 in series with R7)
    R_branch2 = R_p128 + R7
    
    # Calculate the current in the second branch
    I_branch2 = V2 / R_branch2
    
    # Calculate the total current
    I_total = I_branch1 + I_branch2

    print("To find the total current, we can express it as the sum of currents in the two main parallel branches.")
    print("The equation with the given values is:")
    print(f"I_total = ({V2} V / {R3} Ω) + ({V2} V / ( (1 / (1/{R1} Ω + 1/{R2} Ω + 1/{R8} Ω)) + {R7} Ω ))")
    print("\nCalculating the parts of the equation:")
    print(f"Current in Branch 1 = {V2} V / {R3} Ω = {I_branch1:.4f} A")
    print(f"Current in Branch 2 = {V2} V / ({R_p128:.4f} Ω + {R7} Ω) = {I_branch2:.4f} A")
    print("\nSumming the currents:")
    print(f"I_total = {I_branch1:.4f} A + {I_branch2:.4f} A = {I_total:.4f} A")


solve_circuit()