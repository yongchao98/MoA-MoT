import math

# The closed form of the infinite product is derived to be 2^(1/8) * exp(-pi/24),
# which can also be written as (2^(1/8)) / exp(pi/24).

# Let's define the components of this expression.
base_of_power = 2
exponent_numerator = 1
exponent_denominator = 8
pi_constant = math.pi
denominator_of_pi = 24

# The final expression is: base_of_power**(exponent_numerator/exponent_denominator) * math.exp(-pi_constant/denominator_of_pi)

print("The closed-form expression for the infinite product is:")
print(f"({base_of_power}^({exponent_numerator}/{exponent_denominator})) * e^(-{pi_constant:.4f}/{denominator_of_pi})")
print("\nWhich is more simply written as: 2**(1/8) * exp(-pi/24)")

print("\nHere are the numbers in the final equation:")
print(f"Base of the power term: {base_of_power}")
print(f"Exponent of the power term: {exponent_numerator}/{exponent_denominator}")
print(f"The constant pi: ~{pi_constant}")
print(f"The denominator in the exponent of e: {denominator_of_pi}")

# Calculate the numerical value of the expression
numerical_value = base_of_power**(exponent_numerator/exponent_denominator) * math.exp(-pi_constant/denominator_of_pi)

print("\nThe numerical value is:")
print(numerical_value)