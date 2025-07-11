import math

# Define the component values from the circuit diagram
R76 = 76
R8 = 8
R14 = 14
R11 = 11
R29 = 29
V = 41

# --- Calculations ---

# Step 1: Equivalent resistance of the left parallel block (RL)
# R76 and R8 are in parallel
RL = (R76 * R8) / (R76 + R8)

# Step 2: Equivalent resistance of the right block (RR)
# R14 and R11 are in series
RS_right = R14 + R11
# The series combination (RS_right) is in parallel with R29
RR = (RS_right * R29) / (RS_right + R29)

# Step 3: Total equivalent resistance (Req)
# The left block (RL) and right block (RR) are in series
Req = RL + RR

# Step 4: Calculate the total current I
# Using Ohm's Law: I = V / Req
I_final = V / Req

# --- Output ---

print("To find the current I, we first find the total equivalent resistance (Req).")

print("\n1. The two resistors on the left, R_76 and R_8, are in parallel.")
print(f"   RL = ({R76} * {R8}) / ({R76} + {R8}) = {RL:.4f} Ohms.")

print("\n2. On the right, R_14 and R_11 are in series, and this is in parallel with R_29.")
print(f"   Series part (RS) = {R14} + {R11} = {RS_right} Ohms.")
print(f"   RR = ({RS_right} * {R29}) / ({RS_right} + {R29}) = {RR:.4f} Ohms.")

print("\n3. The left and right blocks are in series.")
print(f"   Req = RL + RR = {RL:.4f} + {RR:.4f} = {Req:.4f} Ohms.")

print("\n4. Finally, we use Ohm's Law to find the current I.")
print(f"   The equation is: I = V / Req")
print(f"   Substituting the values: I = {V} V / {Req:.4f} Ohms")
print(f"   I = {I_final:.4f} A")

# Final answer
print(f"\nThe value of the current I is {I_final} A.")