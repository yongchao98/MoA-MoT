import math

# Step 1: Find the limit of the infimum of cos(alpha) as theta -> 0.
# The expression for inf(cos(alpha)) as a function of theta is:
# C(theta) = (4*cos(theta/2)) / sqrt(26 - 10*cos(theta))

# As theta approaches 0:
# cos(theta/2) approaches cos(0) = 1
# cos(theta) approaches cos(0) = 1

# Numerator calculation
num_val_at_0 = 4 * math.cos(0)
print(f"Limit of the numerator 4*cos(theta/2) as theta -> 0 is: {num_val_at_0}")

# Denominator calculation
den_val_at_0 = math.sqrt(26 - 10 * math.cos(0))
print(f"Limit of the denominator sqrt(26 - 10*cos(theta)) as theta -> 0 is: {den_val_at_0}")

# The limit of the infimum of cos(alpha)
limit_cos_alpha = num_val_at_0 / den_val_at_0
print(f"The limit of inf(cos(alpha)) as theta -> 0 is: {limit_cos_alpha}")

# Step 2: Calculate the limit of M(theta).
# M(theta) = arccos(inf(cos(alpha)))
# The limit of M(theta) is arccos of the limit of inf(cos(alpha)).
limit_M_theta = math.acos(limit_cos_alpha)
print(f"The limit of M(theta) is arccos({limit_cos_alpha}) = {limit_M_theta} radians.")

final_answer = limit_M_theta
# The final result is a single number
print(f"The final answer is: {final_answer}")