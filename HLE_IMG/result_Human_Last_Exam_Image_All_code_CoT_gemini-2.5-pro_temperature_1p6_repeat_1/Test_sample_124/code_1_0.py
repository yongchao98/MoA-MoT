import sys

def solve_circuit():
    """
    Calculates the total current flowing through the provided circuit.
    """
    # 1. Define component values from the circuit diagram
    V_source = 1.0  # Voltage in Volts
    R1 = 3.0        # Resistance in Ohms
    R2 = 7.0        # Resistance in Ohms
    R3 = 9.0        # Resistance in Ohms
    R7 = 100.0      # Resistance in Ohms
    R8 = 100.0      # Resistance in Ohms

    # 2. As determined by the analysis, the circuit consists of four parallel branches.
    # We calculate the current through each branch using Ohm's Law (I = V/R).

    # Current through R1
    I1 = V_source / R1

    # Current through R2
    I2 = V_source / R2

    # Current through R3
    I3 = V_source / R3

    # Current through the series combination of R7 and R8
    R78 = R7 + R8
    I78 = V_source / R78

    # 3. The total current is the sum of the currents in the parallel branches.
    I_total = I1 + I2 + I3 + I78

    # 4. Print the final equation with all numbers and the result.
    print("The total current (I_total) is the sum of the currents through the four parallel branches:")
    print("I_total = I_R1 + I_R2 + I_R3 + I_(R7+R8)")
    print(f"I_total = (V / R1) + (V / R2) + (V / R3) + (V / (R7 + R8))")
    # Using sys.stdout.write to print the final equation with each number as requested
    sys.stdout.write("I_total = ({}V / {}Ω) + ({}V / {}Ω) + ({}V / {}Ω) + ({}V / ({}Ω + {}Ω))\n".format(
        V_source, R1, V_source, R2, V_source, R3, V_source, R7, R8))
    sys.stdout.write("I_total = {:.4f} A + {:.4f} A + {:.4f} A + {:.4f} A\n".format(
        I1, I2, I3, I78))

    print(f"\nTotal current = {I_total:.4f} A")
    # Return final numerical answer in the required format
    sys.stdout.write(f"\n<<<{I_total:.4f}>>>\n")

solve_circuit()