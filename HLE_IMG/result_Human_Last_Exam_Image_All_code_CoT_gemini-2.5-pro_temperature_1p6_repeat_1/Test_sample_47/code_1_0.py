import numpy as np

# Define the component values from the circuit diagram
Vp = 41.0
R76 = 76.0
R8 = 8.0
R14 = 14.0
R11 = 11.0
R29 = 29.0
Vg = 0.0

# As determined by the circuit analysis, Node A is at ground potential.
Va = 0.0

print("Step-by-step solution to find the current I:")
print("\n1. Based on the circuit diagram, Node A (the junction of R76, R8, R14) is connected directly to ground.")
print(f"   Therefore, the voltage at Node A is Va = {Va:.1f} V.")

print("\n2. Calculate the voltage at Node B (Vb) using Kirchhoff's Current Law (KCL).")
print("   The KCL equation at Node B is: (Vb - Vp)/R11 + (Vb - Vg)/R29 + (Vb - Va)/R14 = 0")
print(f"   Substituting values: (Vb - {Vp:.0f})/{R11:.0f} + (Vb - {Vg:.0f})/{R29:.0f} + (Vb - {Va:.0f})/{R14:.0f} = 0")
# To solve for Vb, we rearrange the equation: Vb * (1/R11 + 1/R29 + 1/R14) = Vp/R11
Vb = (Vp / R11) / (1/R11 + 1/R29 + 1/R14)
print(f"   Solving for Vb, we find Vb = {Vb:.4f} V.")

print("\n3. Calculate the current I by applying KCL at Node A.")
print("   The current I flows from the ground wire INTO Node A. We first find the sum of all currents entering Node A from the resistors:")
# Current from R76 into A: (Vp - Va) / R76
current_in_76 = (Vp - Va) / R76
# Current from R8 into A: (Vp - Va) / R8
current_in_8 = (Vp - Va) / R8
# Current from R14 into A: (Vb - Va) / R14
current_in_14 = (Vb - Va) / R14

# The current leaving Node A to ground (I_out) must equal the sum of currents entering from resistors.
I_out = current_in_76 + current_in_8 + current_in_14
print(f"   Current from R76 into A = ({Vp:.0f} V - {Va:.0f} V) / {R76:.0f} \u03A9 = {current_in_76:.4f} A")
print(f"   Current from R8 into A  = ({Vp:.0f} V - {Va:.0f} V) / {R8:.0f} \u03A9 = {current_in_8:.4f} A")
print(f"   Current from R14 into A = ({Vb:.4f} V - {Va:.0f} V) / {R14:.0f} \u03A9 = {current_in_14:.4f} A")
print(f"   The total current leaving Node A to ground is {I_out:.4f} A.")

# The current I is defined as flowing INTO Node A from ground, so it's the negative of the current leaving.
I = -I_out
print("\n4. Final Calculation for I.")
print(f"   The arrow for I points into Node A, so I = -(Current leaving Node A).")
print(f"   I = - (Current from R76 + Current from R8 + Current from R14)")
print(f"   I = - ({current_in_76:.4f} A + {current_in_8:.4f} A + {current_in_14:.4f} A)")
print(f"   I = {I:.3f} A")