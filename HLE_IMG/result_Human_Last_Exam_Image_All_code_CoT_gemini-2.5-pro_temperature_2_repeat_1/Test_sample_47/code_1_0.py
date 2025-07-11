import sys

def solve_circuit():
    """
    Solves the circuit problem based on nodal analysis.
    """
    # Define the given values for voltage and resistances
    V_S = 41  # Voltage of the source in Volts
    R76 = 76  # Resistance in Ohms
    R8 = 8    # Resistance in Ohms
    R14 = 14  # Resistance in Ohms
    R29 = 29  # Resistance in Ohms
    V_B = 0   # Voltage at Node B is 0V (grounded)

    # Step 1: Solve for the voltage at Node A (V_A) using KCL.
    # The KCL equation at Node A is: V_A/R76 + (V_A - V_S)/R8 + (V_A - V_B)/R14 = 0
    # Rearranging to solve for V_A: V_A * (1/R76 + 1/R8 + 1/R14) = V_S/R8 + V_B/R14
    # Since V_B = 0, the equation simplifies.
    
    term_vs_r8 = V_S / R8
    coeff_va = (1/R76) + (1/R8) + (1/R14)
    V_A = term_vs_r8 / coeff_va
    
    # Step 2: Calculate the currents flowing from the source.
    # The total current I is the sum of the currents through R8 and R29.
    I_R8 = (V_S - V_A) / R8
    I_R29 = (V_S - V_B) / R29
    
    I_total = I_R8 + I_R29

    # Step 3: Print the final equation with all the numerical values.
    print(f"The equation for the total current I is: I = (V_S - V_A) / R8 + (V_S - V_B) / R29")
    print(f"Intermediate calculation for V_A: {V_A:.2f} V")
    print("Plugging in the values:")
    print(f"I = ({V_S} - {V_A:.2f}) / {R8} + ({V_S} - {V_B}) / {R29}")
    print(f"I = {I_R8:.2f} A + {I_R29:.2f} A")
    print(f"I = {I_total:.2f} A")

solve_circuit()