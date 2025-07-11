import sympy

# Define the variable n for symbolic computation
n = sympy.Symbol('n')

# Define the function f(n) based on the known result from mathematical literature.
# f(n) = 1 + n*(n-1)/2
f_n = 1 + n*(n - 1) / 2

# Define the expression for which we need to find the limit.
# The expression is f(n) / (n * log2(n))
# sympy.log(n, 2) represents log base 2 of n.
expression = f_n / (n * sympy.log(n, 2))

# Calculate the limit of the expression as n approaches infinity
limit_value = sympy.limit(expression, n, sympy.oo)

# The problem asks to output the numbers in the final equation.
# The final equation is: Limit = infinity
# We will print the result obtained from sympy.
print("The limit is:")
print(limit_value)