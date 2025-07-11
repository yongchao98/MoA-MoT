import sympy

# Define symbolic variables
n, eta = sympy.symbols('n eta', positive=True)
p = sympy.Function('p')(n, eta)

# The d-threshold for Hamiltonicity is given by the following asymptotic relation:
# We use log for the natural logarithm as is standard in this field.
threshold_expression = sympy.O(sympy.log(n / eta) / eta)

# Create the equation string
equation_str = f"p(n, eta) = {threshold_expression}"

# The prompt asks to print numbers, but the threshold is an asymptotic formula.
# There is a hidden constant, C, in the Theta notation, so p(n, eta) = C * log(n/eta)/eta.
# Since C is not specified, we present the asymptotic formula.
# Also, there are no other numbers in the final equation derived.
print("The d-threshold for Hamiltonicity is given by the formula:")
print(equation_str)
