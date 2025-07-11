import math

# Step 1: Define the values of k for the two families of solutions.
k1 = 13
k2 = 5

# Step 2: Calculate the characteristic roots for each family.
# The general form is phi_k = (k + sqrt(k^2 - 4)) / 2.
d1 = k1**2 - 4
phi1_val = (k1 + math.sqrt(d1)) / 2

d2 = k2**2 - 4
phi2_val = (k2 + math.sqrt(d2)) / 2

# Step 3: Calculate the limit L = 2/ln(phi_13) + 2/ln(phi_5).
term1 = 2 / math.log(phi1_val)
term2 = 2 / math.log(phi2_val)
limit_L = term1 + term2

# Step 4: Calculate the final result required by the problem.
result = 10**4 * limit_L

# Step 5: Output the details of the calculation and the final answer.
print("The calculation is based on the formula for the limit L:")
print(f"L = 2 / ln(({k1} + sqrt({d1})) / 2) + 2 / ln(({k2} + sqrt({d2})) / 2)")
print("\n--- Intermediate Values ---")
print(f"For k = {k1}:")
print(f"  phi_13 = ({k1} + sqrt({d1})) / 2 = {phi1_val}")
print(f"  Term 1 = 2 / ln(phi_13) = {term1}")
print(f"For k = {k2}:")
print(f"  phi_5 = ({k2} + sqrt({d2})) / 2 = {phi2_val}")
print(f"  Term 2 = 2 / ln(phi_5) = {term2}")
print(f"\nTotal Limit L = {limit_L}")
print("\n--- Final Answer ---")
print(f"The value of 10^4 * L is: {result}")
print(f"The integer part of this value is: {math.floor(result)}")