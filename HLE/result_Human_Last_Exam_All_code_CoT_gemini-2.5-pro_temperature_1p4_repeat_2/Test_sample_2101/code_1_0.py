import math

# Step 1: Define the angles of the isosceles right triangle in radians.
# The vertices are A, B, C. Let B be the vertex with the right angle.
angle_A = math.pi / 4  # 45 degrees
angle_B = math.pi / 2  # 90 degrees
angle_C = math.pi / 4  # 45 degrees

# Step 2: Calculate the expected angle subtended by the hypotenuse.
# The hypotenuse is the side opposite the right angle vertex (B).
# According to Sylvester's formula from geometric probability,
# E[angle subtended by side opposite V] = pi - angle(V).
expected_angle_hypotenuse = math.pi - angle_B

# Step 3: Calculate the probability of escaping through the hypotenuse (P_H).
# This is the expected subtended angle divided by the total angle (2 * pi).
prob_hypotenuse = expected_angle_hypotenuse / (2 * math.pi)

# Step 4: Calculate the final difference (D).
# D = P_H - (P_L1 + P_L2)
# Since P_H + P_L1 + P_L2 = 1, we have P_L1 + P_L2 = 1 - P_H.
# So, D = P_H - (1 - P_H) = 2 * P_H - 1.
difference = 2 * prob_hypotenuse - 1

# Print the final result and the numbers used in the calculation.
print("The angles of the triangle are: pi/2, pi/4, pi/4.")
print(f"The expected angle subtended by the hypotenuse is pi - pi/2 = {expected_angle_hypotenuse/math.pi:.2f}*pi radians.")
print(f"The probability of escaping through the hypotenuse is ({expected_angle_hypotenuse/math.pi:.2f}*pi) / (2*pi) = {prob_hypotenuse}.")
print("The difference between the probability of escaping through the hypotenuse and through the legs is P_H - (1 - P_H) = 2*P_H - 1.")
print(f"Final difference = 2 * {prob_hypotenuse} - 1 = {difference}")

<<< -0.5 >>>