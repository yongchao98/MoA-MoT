import math

# Step 1: Define the active power for each load in MW
load_l1 = 1.85
load_l2 = 1.7
load_l3 = 1.75
load_l4 = 1.9
load_l5 = 2.4

# Step 2: Calculate the total active power (P_total)
loads = [load_l1, load_l2, load_l3, load_l4, load_l5]
p_total = sum(loads)

# Step 3: Define the power factor
power_factor = 0.95

# Step 4: Calculate the required reactive power compensation (Q_compensation)
# The angle phi is calculated from the power factor
phi = math.acos(power_factor)
# The total reactive power is P_total * tan(phi)
q_compensation = p_total * math.tan(phi)

# Step 5: Print the results, showing the equation as requested
print("Step 1: Calculate the total active power (P_total).")
print(f"P_total = {load_l1} MW + {load_l2} MW + {load_l3} MW + {load_l4} MW + {load_l5} MW")
print(f"P_total = {p_total:.2f} MW\n")

print("Step 2: Calculate the required reactive power compensation (Q_compensation).")
print("The required compensation is the total reactive power demand of the loads.")
print(f"The formula is: Q_compensation = P_total * tan(acos(Power_Factor))")
print(f"Using the values:")
print(f"Q_compensation = {p_total:.2f} * tan(acos({power_factor}))")
print(f"Q_compensation = {q_compensation:.2f} MVAR\n")

print("The required reactive power compensation from the ESS at Bus 9 is {:.2f} MVAR.".format(q_compensation))
<<<3.16>>>