from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from A
eq1 = Eq(y - 2, -7/6 * (x + 2))

# Equation of altitude from C
eq3 = Eq(y - 2, -4/3 * (x + 3))

# Solve for intersection of eq1 and eq3
intersection_AC = solve((eq1, eq3), (x, y))

# Check if this intersection satisfies the vertical line equation x = -10
is_on_B_altitude = intersection_AC[x] == -10

# Print the result
print(intersection_AC, is_on_B_altitude)