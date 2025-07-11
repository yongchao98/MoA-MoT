import math
from fractions import Fraction

# Step 1: Define the numerator, which is the binomial coefficient C(50, 25)
# C(n, k) = n! / (k! * (n-k)!)
# For C(50, 25), n=50 and k=25.
numerator = math.comb(50, 25)

# Step 2: Define the denominator, which is 2^50
denominator = 2**50

# Step 3: Create the fraction using the Fraction class
# This keeps the result in its exact fractional form.
result_fraction = Fraction(numerator, denominator)

# The original problem asks for the integral I, which is pi * result_fraction.
# However, the instruction "write the answer strictly as a fraction" implies
# we should provide the rational part of the answer, i.e., I/pi.
# Final Answer format requirement: print the equation
print(f"The integral is equal to (pi * {numerator}) / {denominator}")
print("The value of the integral divided by pi is:")
print(result_fraction)