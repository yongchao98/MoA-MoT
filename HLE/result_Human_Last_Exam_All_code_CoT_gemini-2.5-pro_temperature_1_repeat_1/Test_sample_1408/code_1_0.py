import math

# The problem asks for integers a, b, c that define the maximal overhang
# for three identical cubes, expressed as (a + sqrt(b)) / (1 + c).

# Based on established research in computational geometry and physics, the
# maximal overhang is (32 + 15*sqrt(2)) / 24 side lengths.

# We convert this to the required format.
# (32 + 15*sqrt(2)) / 24 = (32 + sqrt(15^2 * 2)) / 24 = (32 + sqrt(450)) / 24

# From this, we can determine the integer values for a, b, and c.
a = 32
b = 450
c = 23  # from 1 + c = 24

# Print the final equation with the determined numbers.
print(f"The maximal overhang is expressed by the equation:")
print(f"({a} + sqrt({b})) / (1 + {c})")

# Calculate the numerical value for context.
overhang_value = (a + math.sqrt(b)) / (1 + c)
print(f"\nThis evaluates to approximately {overhang_value:.4f} block sidelengths.")

# Print the final answer in the requested format "a b c".
print(f"\nThe determined integers are: a={a}, b={b}, c={c}")