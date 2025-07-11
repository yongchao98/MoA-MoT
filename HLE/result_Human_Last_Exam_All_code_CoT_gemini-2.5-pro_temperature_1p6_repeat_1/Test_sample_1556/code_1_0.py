import math

# Step 3: Define the squared elements of the matrix B based on properties
# a^2 + (a^2)^2 = 1  => a^2 = (sqrt(5)-1)/2
a_squared = (math.sqrt(5) - 1) / 2
b_squared = 1 - a_squared

# For the calculation, we need a and b, not just their squares.
# We assume a and b are positive, without loss of generality.
a = math.sqrt(a_squared)
b = math.sqrt(b_squared)
ab_product = a * b

# Step 4 & 5: Use the experimental data to calculate |α|^2
# From the measured output probabilities P(0)=0.36 and P(1)=0.64,
# the output amplitudes are sqrt(0.36)=0.6 and sqrt(0.64)=0.8.
# alpha_in = 0.6*a + 0.8*b (choosing the positive sign case)
# |alpha_in|^2 = (0.6*a + 0.8*b)^2 = 0.36*a^2 + 0.64*b^2 + 0.96*a*b
term1 = 0.36 * a_squared
term2 = 0.64 * b_squared
term3 = 0.96 * ab_product

alpha_squared = term1 + term2 + term3

print(f"The equation for |α|² is: (0.6 * a + 0.8 * b)²")
print(f"This expands to: 0.36 * a² + 0.64 * b² + 0.96 * a * b")
print(f"Substituting the values:")
print(f"a² = {a_squared}")
print(f"b² = {b_squared}")
print(f"a*b = {ab_product}")
print(f"|α|² = 0.36 * {a_squared} + 0.64 * {b_squared} + 0.96 * {ab_product}")
print(f"|α|² = {term1} + {term2} + {term3}")
print(f"The final value for |α|² is: {alpha_squared}")
