# Step 1: Define the factors for the numerator based on the x-intercepts.
# Root at x=-b and x=b -> (x+b)(x-b) = x^2 - b^2
# Root at x=d (double root) -> (x-d)^2
numerator_part1 = "(x - d)**2"
numerator_part2 = "(x**2 - b**2)"
numerator_expression = f"{numerator_part1} * {numerator_part2}"

# Step 2: Define the factors for the denominator based on the vertical asymptotes.
# Asymptote at x=a (even multiplicity) -> (x-a)^2
# Asymptote at x=c (odd multiplicity) -> (x-c)
denominator_part1 = "(x - a)**2"
denominator_part2 = "(x - c)"
denominator_expression = f"{denominator_part1} * {denominator_part2}"

# Step 3: Combine the parts to form the final equation for f(x).
# The final equation is the ratio of the numerator to the denominator.
print("The equation that represents the graph is:")
print(f"f(x) = ( {numerator_expression} ) / ( {denominator_expression} )")
