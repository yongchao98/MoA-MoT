import math
from fractions import Fraction

# Step 1: Define the weights and number of variables
weights = [22, 29, 49, 50, 75]
n = len(weights)

# Step 2: Calculate the degree 'd' based on the Calabi-Yau condition
d = sum(weights)

# Step 3: Apply the formula for the Crawley-Nordström invariant
# Numerator = (n-2)! * d^(n-3)
n_minus_2_fact = math.factorial(n - 2)
d_pow_n_minus_3 = d ** (n - 3)
numerator = n_minus_2_fact * d_pow_n_minus_3

# Denominator = w1 * w2 * ... * wn
denominator = math.prod(weights)

# Step 4: Calculate the final invariant and simplify the fraction
invariant = Fraction(numerator, denominator)

# Print the calculation steps with the actual numbers
print("The formula for the Crawley-Nordström invariant μ is:")
print("μ = ( (n-2)! * d^(n-3) ) / (w_1 * w_2 * ... * w_n)\n")

print("Substituting the values:")
print(f"n = {n}")
print(f"Weights = {weights}")
print(f"d = sum(weights) = {d}\n")

weight_product_str = " * ".join(map(str, weights))
print(f"μ = ( ({n}-2)! * {d}^({n}-3) ) / ( {weight_product_str} )")
print(f"μ = ( {n_minus_2_fact} * {d_pow_n_minus_3} ) / {denominator}")
print(f"μ = {numerator} / {denominator}\n")

print("The simplified Crawley-Nordström invariant is:")
print(f"μ = {invariant.numerator} / {invariant.denominator}")
