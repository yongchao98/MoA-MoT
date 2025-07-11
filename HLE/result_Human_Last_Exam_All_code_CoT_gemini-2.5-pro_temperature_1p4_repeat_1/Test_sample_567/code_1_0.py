import math

# The question asks for the value of 'a' where the volume constraint becomes the sole
# obstruction for the symplectic embedding of the ellipsoid E(1,a) into a ball.
# According to the work by McDuff and Schlenk, the nature of the obstructions
# changes at a = τ^4, where τ is the golden ratio. For values of a >= τ^4,
# it becomes possible for the embedding capacity c(a) to equal the volume bound sqrt(a).
# This transition value is a = ((1+sqrt(5))/2)^4, which simplifies to (7+3*sqrt(5))/2.

# We define the integer components of the final expression.
numerator_part_1 = 7
numerator_part_2_coeff = 3
number_under_sqrt = 5
denominator = 2

# Calculate the final value of 'a'.
a_value = (numerator_part_1 + numerator_part_2_coeff * math.sqrt(number_under_sqrt)) / denominator

# Print the final equation, showing all the components as requested.
# The result is formatted to a few decimal places for clarity.
print(f"The value 'a' is given by the equation:")
print(f"a = ({numerator_part_1} + {numerator_part_2_coeff} * \u221A{number_under_sqrt}) / {denominator} = {a_value:.15f}")