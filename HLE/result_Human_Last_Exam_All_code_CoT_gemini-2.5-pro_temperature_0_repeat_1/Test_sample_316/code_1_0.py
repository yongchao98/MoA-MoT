import sympy

# Let x = 1/p. The exponent alpha(x) is the maximum of three lines:
# L0(x) = 0
# L1(x) = 1/4 - x  (from transverse/Knapp examples)
# L2(x) = 1/2 - 5x/2 (from planar/Kakeya examples)

# The critical exponents are the p-values where the slope of alpha(x) changes.
# This happens at the intersection points of these lines.

x = sympy.Symbol('x')

# Equation for the first critical point (p=4)
# Intersection of L0 and L1
eq1 = sympy.Eq(0, 1/4 - x)
sol1 = sympy.solve(eq1, x)
p1 = 1/sol1[0]

# Equation for the second critical point
# Intersection of L1 and L2
eq2 = sympy.Eq(1/4 - x, 1/2 - 5*x/2)
sol2 = sympy.solve(eq2, x)
p2 = 1/sol2[0]

# The problem states one critical exponent is 4. We need to find the other one.
print("The exponent alpha(1/p) is determined by the maximum of several linear functions of 1/p.")
print("The relevant lines are:")
print(f"L0(1/p) = 0")
print(f"L1(1/p) = 1/4 - 1/p")
print(f"L2(1/p) = 1/2 - 5/(2*p)")
print("\nThe critical exponents are the p-values where the dominant line changes.")
print(f"The first critical exponent is found by solving L0 = L1:")
print(f"0 = 1/4 - 1/p  => p = {p1}")
print(f"The second critical exponent is found by solving L1 = L2:")
print(f"1/4 - 1/p = 1/2 - 5/(2*p)  => p = {p2}")
print("\nOne critical exponent is given as 4.")
other_exponent = p2 if p1 == 4 else p1
print(f"The other critical exponent is {other_exponent}.")
