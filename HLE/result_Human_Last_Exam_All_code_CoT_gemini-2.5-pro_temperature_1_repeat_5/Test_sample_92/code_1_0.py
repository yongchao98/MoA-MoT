import math

# Step 5: Solve for the roots of 3r^2 - 7r + 3 = 0
# r = (7 +- sqrt(13)) / 6
sqrt13 = math.sqrt(13)
r1 = (7 - sqrt13) / 6
r2 = (7 + sqrt13) / 6

# Step 7: Calculate the final probability
# p0 = (r2^2024 - r1^2024) / (r1 - r2)
# This is equivalent to (r1^2024 - r2^2024) / (r2 - r1)
n = 2024
numerator = r1**n - r2**n
denominator = r2 - r1

probability = numerator / denominator

# The problem is structured as a generalized gambler's ruin problem.
# The derived probability, p_n, represents the probability of reaching state A before state B, starting from n.
# The formula p_0 = (r_1^2024 - r_2^2024) / (r_2 - r_1) correctly represents this probability.
# It can also be expressed as a sum, which is computationally more stable.
# p_0 = sum_{k=0}^{2023} r_2^k * r_1^(2023-k)

# We can compute this sum.
total_sum = 0
power = 2023
for k in range(power + 1):
    term = (r2**k) * (r1**(power-k))
    total_sum += term
    
# Let's print the equation with the values.
# The equation is p0 = (r1^n - r2^n) / (r2 - r1)
print(f"The roots of the characteristic equation are r1 = {r1} and r2 = {r2}.")
print("The probability is calculated using the formula:")
print(f"p(0) = (r1^{n} - r2^{n}) / (r2 - r1)")
print(f"p(0) = (({r1})^{n} - ({r2})^{n}) / (({r2}) - ({r1}))")
print(f"where n = 2024")
print(f"The calculated probability is: {total_sum}")
