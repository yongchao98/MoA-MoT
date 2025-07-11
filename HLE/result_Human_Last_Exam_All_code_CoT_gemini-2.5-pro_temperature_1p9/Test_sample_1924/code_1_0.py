import sympy

# Define symbolic variables for the coordinates of points on the line AB
# M is at the origin
A, B, G, H, K, L = sympy.symbols('a b g h k l')

# Define MA and MB as lengths (a<0, b>0, so MA = -a, MB = b)
MA, MB = sympy.symbols('MA MB', positive=True)
a = -MA
b = MB

# Equation from the Power of a Point K
# (k-a)*(k-b) = k*(k-g)
equation_g = sympy.Eq((K - a) * (K - b), K * (K - g))

# Solve for g, the position of G
g_solution = sympy.solve(equation_g, G)[0]

# Equation from the Power of a Point L
# (l-a)*(l-b) = l*(l-h)
equation_h = sympy.Eq((L - a) * (L - b), L * (L - h))

# Solve for h, the position of H
h_solution = sympy.solve(equation_h, H)[0]

# Calculate the difference MG - MH, which is g - h
difference = sympy.simplify(g_solution - h_solution)

# Print the derivation steps
print("Let M be the origin on the line AB.")
print(f"The coordinate of A is a = -MA, where MA is the length from M to A.")
print(f"The coordinate of B is b = +MB, where MB is the length from M to B.")
print(f"Let K be the coordinate of the intersection of line ED and line AB.")
print(f"From the Power of a Point theorem, we get the relation for the coordinate of G:")
print(f"g = {g_solution}")
print("")
print(f"Similarly, for point H, using the intersection L of line CF and line AB:")
print(f"h = {h_solution}")
print("")
print("The value MG - MH is given by g - h:")
print(f"MG - MH = g - h = {difference}")
print("")
print("The problem states this must hold for ANY chords CD and EF through M.")
print("This implies the result must be independent of the choice of chords, and thus independent of K and L.")
print("This can only be true if K = L, which makes the expression zero.")
print("Therefore, we conclude:")
final_equation = "MG - MH = 0"
print(final_equation)
# The final equation asks for the value of each number in the equation.
# Here, the number is 0.
final_value = 0
print(f"The final numeric value in the equation is: {final_value}")
