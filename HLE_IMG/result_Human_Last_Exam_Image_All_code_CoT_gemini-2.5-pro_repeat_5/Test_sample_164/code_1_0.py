import math

# Step 1: Define the given parameters from the problem description.
# The voltage at Bus 4 is to be maintained at 400 kV.
V_kV = 400.0  # Voltage in kilovolts

# The reactance of the shunt capacitor bank is given as 5 Ohms.
Xc = 5.0      # Reactance in Ohms

# Step 2: Convert the voltage from kilovolts (kV) to volts (V) for the calculation.
V_V = V_kV * 1000

# Step 3: Calculate the reactive power (Q) in VAR using the formula Q = V^2 / Xc.
# V is the voltage in Volts, and Xc is the reactance in Ohms.
Q_VAR = V_V**2 / Xc

# Step 4: Convert the reactive power from VAR to Mega-VAR (MVAR) for standard representation.
# 1 MVAR = 1,000,000 VAR.
Q_MVAR = Q_VAR / 1000000

# Step 5: Print the final answer, clearly showing the equation and the values used.
print("The reactive power (Q_C) supplied by the shunt capacitor bank is calculated using the formula Q_C = V^2 / X_C.")
print(f"To maintain the voltage at Bus 4 at {V_kV} kV, the required reactive power is:")
# The final equation with all numbers is printed below.
print(f"Q_C = ({int(V_V)} V)^2 / {int(Xc)} Ω = {int(Q_VAR)} VAR")
print(f"This is equivalent to {Q_MVAR} MVAR.")

# The final numerical answer.
# print(f"<<<{Q_MVAR}>>>")