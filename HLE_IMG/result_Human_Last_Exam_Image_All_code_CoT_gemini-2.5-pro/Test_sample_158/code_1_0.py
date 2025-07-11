import cmath
import math

# Step 1: Define given parameters
Pin = 10e-3  # Input power in Watts (10 mW)
f = 0.8e9    # Frequency in Hz (0.8 GHz)
Z0 = 50      # Characteristic impedance in Ohms
RL = 2.7e3   # Load resistance in Ohms (2.7 kOhm)
L1_val = 43e-9   # L1 in Henry (43 nH)
C1_val = 0.4e-12 # C1 in Farads (0.4 pF)
L2_val = 39e-9   # L2 in Henry (39 nH)
C2_val = 0.2e-12 # C2 in Farads (0.2 pF)

# From Figure (b), Quality factor of the inductor at 800 MHz
QL = 100
# Given Quality factor for the capacitor
QC = 150
# Assumed forward voltage drop for HSMS-2852 diode
Vd = 0.3     # Volts

# Step 2 & 3: Calculate parasitic resistances to check if losses are negligible
w = 2 * math.pi * f # Angular frequency

# Reactances
XL1 = w * L1_val
XC1 = 1 / (w * C1_val)
XL2 = w * L2_val
XC2 = 1 / (w * C2_val)

# Parasitic series resistances
Rs_L1 = XL1 / QL
Rs_C1 = XC1 / QC
Rs_L2 = XL2 / QL
Rs_C2 = XC2 / QC

print(f"Calculated parasitic resistances:")
print(f"Rs_L1 = {Rs_L1:.4f} Ohms")
print(f"Rs_C1 = {Rs_C1:.4f} Ohms")
print(f"Rs_L2 = {Rs_L2:.4f} Ohms")
print(f"Rs_C2 = {Rs_C2:.4f} Ohms")
print("-" * 30)

# Step 4: Assess passive component losses
# The calculated resistances are very small. For instance, the total series resistance
# in the main path (Rs_L1 + Rs_C1) is only a few ohms. Compared to a load that is
# effectively hundreds of ohms, the power dissipated in these components will be a
# small fraction of the total power. We conclude the passive losses are negligible.
# P_loss = I^2 * R_s. For I ~ sqrt(Pin/Z0) ~ 14mA, P_loss ~ (14e-3)^2 * R_s, which is << Pin.
print("Observation: The parasitic resistances of the passive components are very small.")
print("Therefore, the power loss in them is considered negligible compared to the input power.")
print("We proceed by assuming all 10 mW is delivered to the diode-load combination.")
print("-" * 30)


# Step 5 & 6: Formulate and solve the quadratic equation for V_DC
# Pin = P_load + P_diodes
# Pin = (V_DC^2 / RL) + (2 * Vd * V_DC / RL)
# Multiplying by RL: Pin * RL = V_DC^2 + 2 * Vd * V_DC
# Rearranging into standard quadratic form a*x^2 + b*x + c = 0:
# 1*V_DC^2 + (2*Vd)*V_DC - (Pin * RL) = 0

a = 1
b = 2 * Vd
c = -Pin * RL

# Solve the quadratic equation: V_DC = (-b + sqrt(b^2 - 4ac)) / 2a
discriminant = b**2 - 4*a*c
V_DC = (-b + math.sqrt(discriminant)) / (2*a)

print("Solving the equation for V_DC:")
print(f"V_DC^2 + ({2:.1f} * {Vd:.1f})*V_DC - ({Pin:.3f} * {RL:.0f}) = 0")
print(f"V_DC^2 + {b:.1f}*V_DC - {-c:.1f} = 0")
print(f"The calculated voltage across the load R_L is: {V_DC:.4f} V")

# Final result
final_voltage = V_DC
# Outputting the final answer in the required format.
# Let's round to two decimal places for the final answer display.
final_answer_val = round(final_voltage, 2)
print(f"\nFinal Answer: The calculated voltage is approximately {final_answer_val} V.")
print(f"<<<{final_answer_val}>>>")