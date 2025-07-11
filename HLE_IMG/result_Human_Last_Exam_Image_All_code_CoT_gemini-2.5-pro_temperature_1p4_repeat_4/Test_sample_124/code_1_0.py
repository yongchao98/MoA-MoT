# Define the component values from the circuit diagram
V2 = 1.0  # Voltage of the source in Volts
R1 = 3.0  # Resistance in Ohms
R2 = 7.0  # Resistance in Ohms
R3 = 9.0  # Resistance in Ohms
R7 = 100.0 # Resistance in Ohms
R8 = 100.0 # Resistance in Ohms

# This circuit is a network that requires nodal analysis.
# The wire from the top rail to the junction of R2 and R3 sets that junction's voltage (V_D) to V2.
# The wire between the R1/R2 junction (Node C) and R7/R8 junction (Node E) sets V_C = V_E.
# We solve for the voltage V_C (and V_E) using Kirchhoff's Current Law.

# The equation for V_C is: V_C * (1/R1 + 1/R2 + 1/R7 + 1/R8) = 1/R1 + 1/R2 + 1/R8
numerator = (1/R1) + (1/R2) + (1/R8)
denominator = (1/R1) + (1/R2) + (1/R7) + (1/R8)
V_C = numerator / denominator
V_E = V_C # Since they are connected by a wire

print(f"Calculated voltage at Node C and E: {V_C:.4f} V")

# The total current from the source is the sum of currents returning to ground.
# Current returns to ground via R3 and R7.
V_D = V2 # Due to the wire from the source

# Calculate current through each path to ground
I_R3 = (V_D) / R3
I_R7 = V_E / R7

# The total current is the sum of these two currents.
I_total = I_R3 + I_R7

print("\n--- Calculating Total Current ---")
print(f"The total current is the sum of the currents through R3 and R7.")
print(f"Equation for current through R3: I_R3 = V_D / R3 = {V_D} V / {R3} Ω = {I_R3:.4f} A")
print(f"Equation for current through R7: I_R7 = V_E / R7 = {V_E:.4f} V / {R7} Ω = {I_R7:.4f} A")
print(f"\nFinal Equation for Total Current: I_total = I_R3 + I_R7")
print(f"I_total = {I_R3:.4f} A + {I_R7:.4f} A")
print(f"The total current flowing through the circuit is: {I_total:.4f} A")