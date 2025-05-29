# Equation of altitude from C to AB
eq3 = Eq(y - 7, -11/10 * (x + 3))

# Solve the system of equations
solution2 = solve((eq3, eq2), (x, y))
print(solution2)