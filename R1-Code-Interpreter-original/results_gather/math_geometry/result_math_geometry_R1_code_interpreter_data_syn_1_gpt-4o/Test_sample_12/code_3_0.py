import sympy as sp

# Define the variables
R = 18
r = 2
rs = sp.symbols('rs', positive=True)

# Equation for the slant height
l = R - r + 2*rs

# Pythagorean theorem for the slant height
h = sp.sqrt(l**2 - (R - r)**2)

# Solve for the radius of the sphere
equation = sp.Eq(l, sp.sqrt((R - r)**2 + h**2))
solution = sp.solve(equation, rs)

# Filter the positive solution
positive_solution = [sol.evalf() for sol in solution if sol.is_real and sol > 0]

# Print the positive solution
print(positive_solution[0] if positive_solution else "No solution found")