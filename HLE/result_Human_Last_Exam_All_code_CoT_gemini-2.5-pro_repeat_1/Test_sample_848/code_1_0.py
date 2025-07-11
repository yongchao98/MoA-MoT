import math

# Step 1: Define the constants L_A and L_B derived from analyzing the sequences.
L_A = 13
L_B = 5

# Step 2: Calculate the growth rates phi_A and phi_B.
# These are the larger roots of the characteristic equations z^2 - L*z + 1 = 0.
phi_A = (L_A + math.sqrt(L_A**2 - 4)) / 2
phi_B = (L_B + math.sqrt(L_B**2 - 4)) / 2

# Step 3: Calculate the value of the limit C = lim_{N->inf} F(N)/ln(N).
# The formula is C = 2/ln(phi_A) + 2/ln(phi_B).
limit_C = 2 / math.log(phi_A) + 2 / math.log(phi_B)

# Step 4: Calculate the final result as requested by the problem.
result = 10000 * limit_C

# Output the intermediate values and the final answer.
print(f"L_A = {L_A}")
print(f"L_B = {L_B}")
print(f"phi_A = (13 + sqrt(165))/2 = {phi_A}")
print(f"phi_B = (5 + sqrt(21))/2 = {phi_B}")
print(f"Limit C = 2/ln(phi_A) + 2/ln(phi_B) = {limit_C}")
print(f"Final value = 10^4 * C = {result}")
print(f"The integer part is: {int(result)}")