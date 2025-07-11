import sys

def calculate_circuit_current():
    """
    This function calculates the total current flowing from the voltage source in the given circuit.
    """
    # Define the component values from the circuit diagram
    V2 = 1.0  # Voltage in Volts
    R1 = 3.0  # Resistance in Ohms (bypassed)
    R2 = 7.0  # Resistance in Ohms (bypassed)
    R3 = 9.0  # Resistance in Ohms
    R7 = 100.0 # Resistance in Ohms
    R8 = 100.0 # Resistance in Ohms

    # Step 1: Explain the circuit simplification
    print("Circuit Analysis:")
    print(f"A wire bypasses resistors R1 ({R1} Ω) and R2 ({R2} Ω), so they do not affect the circuit.")
    print("The simplified circuit consists of two parallel branches connected to the 1V source:")
    print(f" - Branch A: Resistor R3 ({R3} Ω)")
    print(f" - Branch B: Resistors R7 ({R7} Ω) and R8 ({R8} Ω) in series.")
    print("-" * 40)

    # Step 2: Calculate the equivalent resistance of the series branch (R7 and R8)
    R_series_78 = R7 + R8
    print("Calculation Steps:")
    print(f"1. Resistance of the series branch (R7 + R8):")
    print(f"   R_series = {R7} Ω + {R8} Ω = {R_series_78} Ω")
    
    # Step 3: Calculate the total equivalent resistance of the parallel circuit
    R_total = (R3 * R_series_78) / (R3 + R_series_78)
    print("\n2. Total equivalent resistance (R3 in parallel with R_series):")
    print(f"   R_total = ({R3} Ω * {R_series_78} Ω) / ({R3} Ω + {R_series_78} Ω) = {R_total:.4f} Ω")

    # Step 4: Calculate the total current using Ohm's Law
    I_total = V2 / R_total
    print("\n3. Total current from the source (Ohm's Law: I = V / R):")
    print(f"   I_total = {V2} V / {R_total:.4f} Ω")

    print("-" * 40)
    print(f"The total current flowing through the circuit is {I_total:.4f} A.")

if __name__ == "__main__":
    calculate_circuit_current()
