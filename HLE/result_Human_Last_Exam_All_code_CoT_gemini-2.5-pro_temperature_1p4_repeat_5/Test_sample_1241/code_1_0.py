import numpy as np

# Step 1: Define the given transition rates
l01 = 0.019
l10 = 0.65
l12 = 0.4
l21 = 0.392
l23 = 0.008
l31 = 0.008

# Step 2: Set up the system of linear equations for the steady state (p_i' = 0)
# We can express p0, p2, and p3 in terms of p1.
# From 0 = -l01*p0 + l10*p1  => p0 = (l10/l01)*p1
# From 0 = l12*p1 - (l21+l23)*p2 => p2 = (l12/(l21+l23))*p1
# From 0 = l23*p2 - l31*p3 => p3 = (l23/l31)*p2

# Calculate the coefficients for p0, p2, p3 in terms of p1
coeff_p0 = l10 / l01
coeff_p2 = l12 / (l21 + l23)
coeff_p3 = (l23 / l31) * coeff_p2

# Step 3: Use the normalization condition p0 + p1 + p2 + p3 = 1 to solve for p1
# p1*(coeff_p0 + 1 + coeff_p2 + coeff_p3) = 1
sum_of_coeffs = coeff_p0 + 1 + coeff_p2 + coeff_p3
p1 = 1 / sum_of_coeffs

# Step 4: Calculate p0 using the relation from Step 2
p0 = coeff_p0 * p1

# Step 5: Calculate the final required sum
result = p0 + p1

# Output the results in the required format
print(f"By solving the system of linear equations for the steady-state probabilities, we find:")
print(f"P_0(+infinity) = {p0}")
print(f"P_1(+infinity) = {p1}")
print(f"\nThus, the value of P_0(+infinity) + P_1(+infinity) is:")
print(f"{p0} + {p1} = {result}")