# This script formats and prints the final symbolic formula for the given infinite product.

# The left-hand side of the equation, as given in the problem.
# The numbers 3 and 3 are present here.
lhs = "Product_{n=3 to infinity} (1 - z^3/n^3)"

# The denominator terms that arise from splitting the product.
# They correspond to n=1 and n=2. We write 1^3 and 2^3 to show their origin.
term_n1 = "(1 - z^3/(1^3))"
term_n2 = "(1 - z^3/(2^3))"

# The product of Gamma functions derived from the infinite product starting at n=1.
# The numbers 1, 2, 3, 4 are present in this part of the expression.
gamma_product = "Gamma(1 - z) * Gamma(1 - z*exp(i*2*pi/3)) * Gamma(1 - z*exp(i*4*pi/3))"

# Construct the right-hand side of the equation.
rhs = f"1 / ({term_n1} * {term_n2} * {gamma_product})"

# Print the final equation, showing how all the numbers come together.
print("The final expression for the infinite product is:")
print(f"{lhs} = {rhs}")