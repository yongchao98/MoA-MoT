import math

# Define the given values for loads and power factor
l1_mw = 1.85
l2_mw = 1.7
l3_mw = 1.75
l4_mw = 1.9
l5_mw = 2.4
power_factor = 0.95

# Step 1: Calculate the total real power (P_total) by summing all loads
p_total = l1_mw + l2_mw + l3_mw + l4_mw + l5_mw

# Step 2: Calculate the total reactive power (Q_total)
# The angle phi is the arccosine of the power factor
phi = math.acos(power_factor)

# The reactive power Q is P * tan(phi)
q_total = p_total * math.tan(phi)

# Step 3: Print the results and the final equation
# The required reactive power from the ESS is equal to the total reactive power demand.
print("Step 1: Calculate the total real power demand (P_total).")
print(f"P_total = {l1_mw} MW + {l2_mw} MW + {l3_mw} MW + {l4_mw} MW + {l5_mw} MW")
print(f"P_total = {p_total:.2f} MW\n")

print("Step 2: Calculate the total reactive power demand (Q_total) using the power factor.")
print(f"The equation for reactive power is: Q_total = P_total * tan(acos(PF))")
print(f"Q_total = {p_total:.2f} * tan(acos({power_factor}))")
print(f"Q_total = {q_total:.2f} Mvar\n")

print("The required reactive power compensation from the ESS is {:.2f} Mvar.".format(q_total))
print("<<<3.16>>>")