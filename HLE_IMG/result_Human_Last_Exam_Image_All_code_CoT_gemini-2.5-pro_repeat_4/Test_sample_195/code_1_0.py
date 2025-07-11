import sympy

# Define the variable x and the constants a, b, c, d
x, a, b, c, d = sympy.symbols('x a b c d')

# Numerator based on the roots of the function
# Root at -b: (x+b)
# Root at +b: (x-b)
# Double root at d: (x-d)^2
numerator = (x**2 - b**2) * (x - d)**2

# Denominator based on the vertical asymptotes
# Asymptote at a (odd power): (x-a)
# Asymptote at c (even power): (x-c)^2
denominator = (x - a) * (x - c)**2

# The function f(x)
f_x = numerator / denominator

# Print the equation for f(x)
# The prompt asks to "output each number in the final equation". As we are using symbolic constants,
# we will print the equation in its symbolic form.
equation_string = f"f(x) = ({sympy.pretty(numerator, use_unicode=False)}) / ({sympy.pretty(denominator, use_unicode=False)})"

# To make the output cleaner and more like standard math notation, we will format it.
# The following code prints the final equation in a readable format.
# We print the components to satisfy the "output each number" constraint symbolically.

print("The equation for the function f(x) is:")
print(f"f(x) = ( (x**2 - b**2) * (x - d)**2 ) / ( (x - a) * (x - c)**2 )")
print("\nWhere the terms are constructed from the graph's features:")
print(f"Numerator term for roots at -b and b: (x**2 - b**2)")
print(f"Numerator term for the double root at d: (x - d)**2")
print(f"Denominator term for the asymptote at a: (x - a)")
print(f"Denominator term for the asymptote at c: (x - c)**2")
# The values 1 and 2 are the powers.
print(f"The implied powers are 1, 2, 1, 2.")
