# Define the given values from the circuit diagram
V_s = 41.0  # Voltage of the source in Volts
R76 = 76.0  # Resistance in Ohms
R8 = 8.0    # Resistance in Ohms
R11 = 11.0  # Resistance in Ohms
R29 = 29.0  # Resistance in Ohms
R14 = 14.0  # Bridge resistance in Ohms

# --- Step 1: Calculate the voltage Va at the junction of R76 and R8 ---
# This is a voltage divider. Va = V_s * R8 / (R76 + R8)
Va = V_s * R8 / (R76 + R8)

# --- Step 2: Calculate the voltage Vb at the junction of R11 and R29 ---
# This is another voltage divider. Vb = V_s * R29 / (R11 + R29)
Vb = V_s * R29 / (R11 + R29)

# --- Step 3: Calculate the current I flowing through R14 ---
# I = (Va - Vb) / R14
I = (Va - Vb) / R14

# --- Print the results step-by-step ---
print("This circuit is analyzed as a Wheatstone bridge.")
print("The current I is the current flowing through the 14 Ohm resistor.")
print("\nStep 1: Calculate the voltage Va at the midpoint of the left branch.")
print(f"Va = V_s * R8 / (R76 + R8)")
print(f"Va = {V_s} V * {R8} \u03A9 / ({R76} \u03A9 + {R8} \u03A9)")
print(f"Va = {V_s * R8} / {R76 + R8} = {Va:.4f} V")

print("\nStep 2: Calculate the voltage Vb at the midpoint of the right branch.")
print(f"Vb = V_s * R29 / (R11 + R29)")
print(f"Vb = {V_s} V * {R29} \u03A9 / ({R11} \u03A9 + {R29} \u03A9)")
print(f"Vb = {V_s * R29} / {R11 + R29} = {Vb:.4f} V")

print("\nStep 3: Calculate the current I.")
print(f"I = (Va - Vb) / R14")
print(f"I = ({Va:.4f} V - {Vb:.4f} V) / {R14} \u03A9")
print(f"I = {Va - Vb:.4f} V / {R14} \u03A9")
print(f"\nThe value of the current I is: {I:.4f} A")
print("\nThe negative sign indicates that the current flows from right to left, opposite to the direction of the arrow in the diagram.")