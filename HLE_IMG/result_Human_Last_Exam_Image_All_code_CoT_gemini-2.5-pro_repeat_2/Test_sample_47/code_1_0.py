import math

def solve_circuit():
    """
    Calculates the current I in the given circuit diagram.
    """
    # Define the component values
    V = 41.0   # Voltage source in Volts
    R76 = 76.0 # Resistor values in Ohms
    R8 = 8.0
    R14 = 14.0
    R11 = 11.0
    R29 = 29.0

    print("Step-by-step calculation to find the current I:\n")

    # Step 1: Calculate the equivalent resistance of the top parallel section
    r_top_val = (R76 * R14) / (R76 + R14)
    print(f"1. The top section has two resistors in parallel (76 Ohm and 14 Ohm).")
    print(f"   R_top = ({R76} * {R14}) / ({R76} + {R14}) = {r_top_val:.4f} Ohms\n")

    # Step 2: Calculate the equivalent resistance of the middle parallel section
    r_mid_val = (R8 * R11) / (R8 + R11)
    print(f"2. The middle section has two resistors in parallel (8 Ohm and 11 Ohm).")
    print(f"   R_mid = ({R8} * {R11}) / ({R8} + {R11}) = {r_mid_val:.4f} Ohms\n")
    
    # Step 3: The bottom section is a single resistor
    r_bot_val = R29
    print(f"3. The bottom section is a single resistor of {R29} Ohms.\n")

    # Step 4: Calculate the total equivalent resistance for the series circuit
    r_total_val = r_top_val + r_mid_val + r_bot_val
    print("4. The total resistance is the sum of these three sections in series.")
    print(f"   R_total = R_top + R_mid + R_bot = {r_top_val:.4f} + {r_mid_val:.4f} + {r_bot_val} = {r_total_val:.4f} Ohms\n")

    # Step 5: Calculate the current I using Ohm's Law
    current_I = V / r_total_val
    print("5. Finally, calculate the current I using Ohm's Law (I = V / R_total).")
    print("   The final equation with all the numbers is:")
    print(f"   I = {V} / ( (({R76} * {R14}) / ({R76} + {R14})) + (({R8} * {R11}) / ({R8} + {R11})) + {R29} )")
    print(f"   I = {V} / {r_total_val:.4f}")
    print(f"   I = {current_I:.4f} A\n")

    print(f"The final value of the current I is {current_I:.2f} A.")

solve_circuit()