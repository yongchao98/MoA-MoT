import math

# Step 1: Define constants and values for the terms
pi = math.pi
# The value of the integral is a sum of four terms derived from the analytical solution.
# I = I_1 + I_2 + I_3 + I_4
# I_1 = integral of (p^7)/(e^p-1) from 0 to inf
# I_2 = integral of (p)/(e^p-1) from 0 to inf
# I_3 = integral of (p*e^-p)/(e^p-1) from 0 to inf
# I_4 = integral of (e^(p/4)-e^(-p/4))/(2*(e^p-1)) from 0 to inf

# Step 2: Calculate the symbolic and numeric value of each term.
term1_coeff_num = 8
term1_coeff_den = 15
term1_val = (term1_coeff_num / term1_coeff_den) * (pi**8)

term2_coeff_num = 1
term2_coeff_den = 3
term2_val = (term2_coeff_num / term2_coeff_den) * (pi**2)

term3_coeff_num = -1
term3_coeff_den = 2
term3_val = (term3_coeff_num / term3_coeff_den) * pi

term4_val = 1

final_value = term1_val + term2_val + term3_val + term4_val

# Step 3: Print the equation with all the numbers.
# The user wants each number in the final equation to be outputted.
print("The final result is calculated from the sum of four evaluated integrals:")
print(f"Term 1 (from p^7): (8/15) * pi^8 = {term1_val}")
print(f"Term 2 & part of Term 3 (from p): (1/3) * pi^2 = {term2_val}")
print(f"Term 4 (from sinh): (-1/2) * pi = {term3_val}")
print(f"Constant terms: -1 + 2 = {term4_val}")
print("\nFinal Equation:")
print(f"{term1_coeff_num}/{term1_coeff_den} * pi^8 + {term2_coeff_num}/{term2_coeff_den} * pi^2 - {abs(term3_coeff_num)}/{term3_coeff_den} * pi + {term4_val}")
print("\nFinal Value:")
print(final_value)
