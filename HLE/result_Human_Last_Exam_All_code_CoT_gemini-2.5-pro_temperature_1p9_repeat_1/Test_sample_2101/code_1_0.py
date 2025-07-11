import math

# Step 1: Define the internal angles of the isosceles right triangle in radians.
angle_opposite_hypotenuse = math.pi / 2  # 90 degrees
angle_opposite_leg = math.pi / 4       # 45 degrees

# Step 2: Use the formula P_side = (pi - angle_opposite) / (2 * pi) to find probabilities.

# Probability of escaping through the hypotenuse
prob_hypotenuse = (math.pi - angle_opposite_hypotenuse) / (2 * math.pi)

# Probability of escaping through one of the legs
prob_leg = (math.pi - angle_opposite_leg) / (2 * math.pi)

# Step 3: Calculate the difference requested by the problem.
# The total probability of escaping through the legs is 2 * prob_leg.
difference = prob_hypotenuse - (prob_leg + prob_leg)

# Step 4: Display the calculation and the final answer.
# We present the result as fractions to show the exact values.
# P_H = (pi/2) / (2*pi) = 1/4
# P_L = (3*pi/4) / (2*pi) = 3/8
ph_num, ph_den = 1, 4
pl_num, pl_den = 3, 8

print("The problem is to find: P(hypotenuse) - (P(leg1) + P(leg2))")
print(f"The calculation using the derived probabilities is:")
print(f"{ph_num}/{ph_den} - ({pl_num}/{pl_den} + {pl_num}/{pl_den}) = {difference}")
