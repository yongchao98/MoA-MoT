import sympy as sp

# Define the variables
R = 18
r = 2
rs = sp.symbols('rs', positive=True)

# Equation for the radius of the sphere
# The height h is expressed in terms of the radii and the radius of the sphere
h = sp.sqrt((R - r)**2 + (R + r - 2*rs)**2)

# Equation for the sphere inscribed in the truncated cone
equation = sp.Eq(rs, (R * r * h) / (sp.sqrt((R - r)**2 + h**2) + R + r))

# Solve the equation for rs
solution = sp.solve(equation, rs)

# Filter the positive solution
positive_solution = [sol.evalf() for sol in solution if sol.is_real and sol > 0]

# Print the positive solution
print(positive_solution[0])