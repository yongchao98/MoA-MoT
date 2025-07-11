import sympy

# Define the symbols for time t and the function h(t)
t = sympy.symbols('t')
h = sympy.Function('h')(t)

# Determine the explicit form of H(t) based on the derivation
# H(t) = exp(h(t)/2)
# The numbers in the equation are 1 and 2.
H = sympy.exp(h / 2)

# Create a pretty string representation of the equation for H(t)
# We output each number in the final equation as requested.
# The number 1 is the implicit coefficient of h.
# The number 2 is the denominator.
equation_str = f"H(t) = exp((1 * h(t)) / 2)"

print("The explicit form of H(t) is:")
# We can print the Sympy object for a nice mathematical display
# or print the string we constructed. Both convey the result.
# sympy.init_printing(use_unicode=True)
# print(H)
# For a clear text output:
print(equation_str)

# Extracting the components of the expression to fulfill the "output each number" requirement literally
one = sympy.Integer(1)
two = sympy.Integer(2)
print(f"\nThe numbers in the final equation H(t) = exp(h(t)/2) are {one} and {two}.")