import math

# Step 1: Define the given constants from the problem
AD = 12 # height
AE = 13 # angle bisector length

# Step 2: Calculate the distance DE
# In the right-angled triangle ADE, DE = sqrt(AE^2 - AD^2)
DE = math.sqrt(AE**2 - AD**2)

# Step 3: Establish the relationship between the lengths on the base BC.
# Let f be the distance DF from the foot of the altitude D to the foot of the median F.
# Let P be the product of the lengths of the segments BD and CD, i.e., P = BD * CD.
# The relationship derived from the angle bisector property is P = (119/5) * f - 144.
# (This comes from 119*(CD-BD) = 10*(144 + BD*CD) and f = (CD-BD)/2)

# Step 4: Apply the conditions on P to find the range of f.
# Condition 1: For ∠A to be acute, we need P < 144.
# (119/5) * f - 144 < 144
# (119/5) * f < 288
# f < 288 * 5 / 119
f_upper_bound = 288 * 5 / 119

# Condition 2: For a valid triangle with positive side lengths, we need P > 0.
# (119/5) * f - 144 > 0
# (119/5) * f > 144
# f > 144 * 5 / 119
f_lower_bound = 144 * 5 / 119

# Step 5: Calculate the range of m using the relation m^2 = AD^2 + f^2
# Lower bound for m
m_sq_lower = AD**2 + f_lower_bound**2
m_lower_bound = math.sqrt(m_sq_lower)

# Upper bound for m
m_sq_upper = AD**2 + f_upper_bound**2
m_upper_bound = math.sqrt(m_sq_upper)

# Step 6: Print the results.
# The numbers in the final inequality defining the range of m are m_lower_bound and m_upper_bound.
print("The range for the median m is given by the inequality:")
print(f"{m_lower_bound:.4f} < m < {m_upper_bound:.4f}")
print("\nOr, in exact fractional form:")
# For the upper bound, we use the Pythagorean triple (119, 120, 169)
# m_upper_bound = sqrt(12^2 + (1440/119)^2) = sqrt(12^2 * (1 + (120/119)^2)) = 12 * sqrt((119^2 + 120^2)/119^2) = 12 * 169/119
m_upper_exact_num = 12 * 169
m_upper_exact_den = 119

# For the lower bound, m_lower_bound = sqrt(12^2 + (720/119)^2) = 12 * sqrt(1 + (60/119)^2) = (12/119) * sqrt(119^2 + 60^2)
# 119^2 + 60^2 = 14161 + 3600 = 17761
m_lower_exact_num_sqrt = 17761

print(f"(12 * sqrt({m_lower_exact_num_sqrt})) / {m_upper_exact_den} < m < {m_upper_exact_num} / {m_upper_exact_den}")
