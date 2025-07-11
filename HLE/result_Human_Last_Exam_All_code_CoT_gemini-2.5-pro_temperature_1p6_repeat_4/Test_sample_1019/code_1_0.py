import sympy

# Define symbols
n, eta = sympy.symbols('n eta')

# The size of the two components in the worst-case graph construction
# For simplicity in the symbolic expression, we ignore the ceiling/floor functions,
# as they do not affect the leading terms.
c1_size = n/2 - eta
c2_size = n - 1 - c1_size

# The threshold p is inversely proportional to the product of the component sizes
denominator = c1_size * c2_size

# Expand and simplify the denominator
simplified_denominator = sympy.expand(denominator)

# Create the equation string for p
# We output the expanded form for clarity.
# Note: The problem asks for the d-threshold, which is a probability p.
# Our analysis shows p should be Theta(1/denominator).
# We represent the threshold as 1/denominator.
equation = f"p = 1 / ({simplified_denominator})"

# To satisfy the prompt "output each number in the final equation",
# we can format the printout this way. The numbers are 1, 4, 2, 1, 2.
# Using sympy's output can be complex, so we construct the string from the expanded form.
# simplified_denominator = n**2/4 - n/2 - eta**2 + eta
print("The d-threshold for Hamiltonicity is given by p where:")
print(f"p = 1 / ( (1/4)*n**2 - (1/2)*n - eta**2 + eta )")
