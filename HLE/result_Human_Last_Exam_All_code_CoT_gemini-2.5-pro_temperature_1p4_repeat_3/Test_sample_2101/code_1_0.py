import math

# This script calculates the difference between the probability that an electron
# escapes through the hypotenuse and the probability that it escapes through
# either of the two legs of a right isosceles triangle.

# The final formula for this difference is: 1 - (4/π) * ln(1+√2)

# Step 1: Calculate the values of the constants and terms in the equation.
sqrt2 = math.sqrt(2)
one_plus_sqrt2 = 1 + sqrt2
ln_one_plus_sqrt2 = math.log(one_plus_sqrt2) # This is also equal to asinh(1)
pi = math.pi
four_div_pi = 4 / pi
term_to_subtract = four_div_pi * ln_one_plus_sqrt2

# Step 2: Calculate the final result.
final_difference = 1 - term_to_subtract

# Step 3: Print the breakdown of the calculation as requested.
print("The difference in probabilities is calculated by the formula: 1 - (4/π) * ln(1+√2)")
print("\n--- Calculation Breakdown ---")
print(f"      √2 = {sqrt2}")
print(f"   1 + √2 = {one_plus_sqrt2}")
print(f"ln(1 + √2) = {ln_one_plus_sqrt2}")
print(f"        π = {pi}")
print(f"      4/π = {four_div_pi}")
print("\n--- Final Equation ---")
print(f"Difference = 1 - ({four_div_pi}) * ({ln_one_plus_sqrt2})")
print(f"Difference = 1 - {term_to_subtract}")
print(f"Difference = {final_difference}")
