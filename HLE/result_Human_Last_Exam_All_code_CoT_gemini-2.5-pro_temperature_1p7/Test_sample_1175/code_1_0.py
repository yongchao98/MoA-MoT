import math

# These are the constants from the final simplified expression:
# (3/2) * 10^(10/3) + 37/4
a = 3/2
base = 10
exponent_num = 10
exponent_den = 3
c = 37/4

exponent = exponent_num / exponent_den

# Calculate the final result based on the simplified expression
result = a * (base ** exponent) + c

# As requested, outputting each number in the final simplified equation
print("The final simplified equation has the form: a * base**(num/den) + c")
print(f"The number for 'a' is: {a}")
print(f"The number for 'base' is: {base}")
print(f"The number for the exponent numerator 'num' is: {exponent_num}")
print(f"The number for the exponent denominator 'den' is: {exponent_den}")
print(f"The number for 'c' is: {c}")

# Display the simplified equation with its numerical values
print(f"\nFinal Equation: {a} * {base}**({exponent_num}/{exponent_den}) + {c}")

# Print the final numerical answer
print(f"\nCalculated Value: {result}")