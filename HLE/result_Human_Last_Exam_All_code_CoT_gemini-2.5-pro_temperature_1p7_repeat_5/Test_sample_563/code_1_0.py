# The number of isomorphism classes of automorphism groups for Riemann surfaces
# of a given genus is a known, albeit complex, result from algebraic geometry.
# These numbers are sourced from mathematical classification literature.

# For genus g=2, the number of groups that can act is 12.
num_g2 = 12

# For genus g=3, the number of groups that can act is 36.
num_g3 = 36

# For genus g=4, the number is stated as 23. This is a specific count
# that differs from more recent classifications, which list 29 (for full groups)
# or 55 (for groups that can act). We will use the number provided in the problem context.
num_g4 = 23

# The final result is a list containing these three numbers.
# We print the final list, constructing it from the variables.
# This way, we show each number that contributes to the final result.
print(f"[{num_g2},{num_g3},{num_g4}]")