import math

def calculate_circuit_current():
    """
    Calculates the total current flowing from the voltage source in the given circuit.

    The circuit analysis reveals four parallel branches connected across the voltage source:
    1. A branch with resistor R1.
    2. A branch with resistor R2.
    3. A branch with resistor R3.
    4. A branch with resistors R7 and R8 in series.

    The total current is the sum of the currents through each of these parallel branches.
    """
    # Define the component values from the circuit diagram
    V2 = 1.0  # Voltage in Volts
    R1 = 3.0  # Resistance in Ohms
    R2 = 7.0  # Resistance in Ohms
    R3 = 9.0  # Resistance in Ohms
    R7 = 100.0 # Resistance in Ohms
    R8 = 100.0 # Resistance in Ohms

    # Calculate the current through each parallel branch using Ohm's Law (I = V/R)
    I1 = V2 / R1
    I2 = V2 / R2
    I3 = V2 / R3
    
    # For the fourth branch, the total resistance is the sum of R7 and R8
    R_78 = R7 + R8
    I_78 = V2 / R_78

    # The total current is the sum of the currents from all parallel branches
    I_total = I1 + I2 + I3 + I_78

    # Print the explanation and the final equation with all values
    print("The total current (I_total) is the sum of the currents in the parallel branches.")
    print(f"I_total = V2/R1 + V2/R2 + V2/R3 + V2/(R7 + R8)")
    print(f"I_total = {V2}/{R1} + {V2}/{R2} + {V2}/{R3} + {V2}/({R7} + {R8})")
    print(f"I_total = {I1:.4f} A + {I2:.4f} A + {I3:.4f} A + {I_78:.4f} A")
    
    # Print the final result
    print(f"\nThe total current flowing through the circuit is: {I_total:.4f} A")

# Execute the calculation
calculate_circuit_current()