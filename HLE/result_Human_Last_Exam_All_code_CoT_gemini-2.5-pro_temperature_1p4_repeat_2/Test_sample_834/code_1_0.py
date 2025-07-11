import math

# Step 1: Define the geometric constant
# a = 12^(1/4), so a^4 = 12
a_fourth = 12

# Step 2: Calculate area moments of inertia I_ss and I_zz
# I_zz = 55 * a^4 / 12
I_zz = 55 * a_fourth / 12
# I_ss = 73 * a^4 / 12
I_ss = 73 * a_fourth / 12

print(f"Calculated I_zz = {I_zz}")
print(f"Calculated I_ss = {I_ss}")
print("-" * 20)

# Step 3: Calculate L and q0 using the given data and calculated inertia values
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

print(f"Given L = 30 * I_zz / 73, so L = {L}")
print(f"Given q0 = 9 * I_ss / 55, so q0 = {q0}")
print("-" * 20)

# Step 4: Calculate the required force F using the derived relation F = q0*L/10
# This is the final equation with the calculated numbers
print(f"The final equation for F is: F = q0 * L / 10")
print(f"Substituting the values: F = {q0} * {L} / 10")

F = (q0 * L) / 10

# Step 5: Print the final result
print("-" * 20)
print(f"The required force F is: {F}")

# Final Answer in the requested format
print(f"<<<{F}>>>")