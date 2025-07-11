# Define the given values for voltage and resistances
V_source = 41.0  # V
R76 = 76.0       # Ohm
R8 = 8.0         # Ohm
R14 = 14.0       # Ohm
R11 = 11.0       # Ohm
R29 = 29.0       # Ohm

# Based on the circuit analysis:
Va = 0.0   # Potential at Node A (junction of R76, R8)
Vb = 41.0  # Potential at Node B (junction of R8, R14)

# Step 1: Solve for Vc using KCL at Node C.
# Equation: (Vb - Vc) / R14 = Vc / R11 + Vc / R29
# Rearranging for Vc: Vc = (Vb / R14) / (1/R14 + 1/R11 + 1/R29)
print("Step 1: Calculate the potential at Node C (Vc).")
print(f"The KCL equation at Node C is: ({Vb} - Vc) / {R14} = Vc / {R11} + Vc / {R29}")
G14 = 1/R14
G11 = 1/R11
G29 = 1/R29
Vc = (Vb * G14) / (G14 + G11 + G29)
print(f"The potential at Node C (Vc) is {Vc:.4f} V.\n")

# Step 2: Calculate the current I.
# I is the sum of currents leaving Node B.
# I = I_R8 + I_R14
print("Step 2: Calculate the current I.")
# Current through R8
I_R8 = (Vb - Va) / R8
# Current through R14
I_R14 = (Vb - Vc) / R14
# Total current I
I = I_R8 + I_R14

print(f"The equation for the current I is: I = (Vb - Va) / R8 + (Vb - Vc) / R14")
print(f"I = ({Vb} - {Va}) / {R8} + ({Vb} - {Vc:.4f}) / {R14}")
print(f"I = {I_R8:.4f} A + {I_R14:.4f} A")
print(f"The value of the current I is {I:.4f} A.")
