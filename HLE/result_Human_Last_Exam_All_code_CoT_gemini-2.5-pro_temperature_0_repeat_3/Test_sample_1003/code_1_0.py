import math

# This script calculates the required ratio based on the principles of special relativity.

# Step 1: Define the known quantities from the problem statement.
# In the first frame, for any pair of stars i, j, cos(theta_ij) = -1/3.
one_minus_cos_theta = 4/3

# In the second frame, we are given the angles, from which we find 1 - cos(theta').
one_minus_cos_theta_prime_12 = 1 - 0
one_minus_cos_theta_prime_13 = 1 - (-1 / math.sqrt(2))

# Step 2: Use the relativistic aberration formula to find the terms needed for the final ratio.
# The final ratio simplifies to Num / Den, where:
# Den = gamma * (1 - v.n_1)
# Num = gamma * (1 - v.n_3)

# From the angle between S1 and S2, we find the denominator term 'Den'.
# 1 - cos(theta'_12) = (4/3) / Den^2
# 1 = (4/3) / Den^2  => Den^2 = 4/3
Den = math.sqrt(one_minus_cos_theta / one_minus_cos_theta_prime_12)

# From the angle between S1 and S3, we find the numerator term 'Num'.
# 1 - cos(theta'_13) = (4/3) / (Den * Num)
Num = one_minus_cos_theta / (Den * one_minus_cos_theta_prime_13)

# Step 3: Calculate the final ratio, which is Num / Den.
final_ratio = Num / Den

# Step 4: Print the results, showing the numbers that form the final equation.
print("The problem is to calculate the ratio (1 - cos(theta_14)) / (1 - cos(theta_34)).")
print("This ratio simplifies to a division of two terms, which we call Num and Den.")
print("The final equation is: Ratio = Num / Den")
print("")
print(f"The value of the denominator term 'Den' is: {Den}")
print(f"The value of the numerator term 'Num' is: {Num}")
print("")
print("The final equation with the calculated numbers is:")
print(f"{Num} / {Den} = {final_ratio}")
print("")
# The exact symbolic value is 2 - sqrt(2).
print(f"The exact value is 2 - sqrt(2), which is approximately {2 - math.sqrt(2)}.")
