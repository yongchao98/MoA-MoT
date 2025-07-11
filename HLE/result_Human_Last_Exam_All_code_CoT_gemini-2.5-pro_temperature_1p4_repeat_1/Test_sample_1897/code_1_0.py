import math

# The closed-form expression for the infinite product is 2^(1/8) * exp(-pi/24).
# This script will print the components of this expression and its numerical value.

base = 2
exponent_numerator = 1
exponent_denominator = 8
exp_numerator = -1
exp_denominator = 24

print("The closed form of the infinite product is an equation with the following numbers:")
print(f"Base: {base}")
print(f"Exponent: {exponent_numerator}/{exponent_denominator}")
print(f"Base of the exponential: e")
print(f"Exponent of e: pi * ({exp_numerator}/{exp_denominator})")

print("\nThe final expression is:")
print(f"{base}**({exponent_numerator}/{exponent_denominator}) * exp(pi * ({exp_numerator}/{exp_denominator}))")

# Calculate the numerical value
value = math.pow(base, exponent_numerator / exponent_denominator) * math.exp(math.pi * exp_numerator / exp_denominator)

print("\nWhich evaluates to:")
print(value)