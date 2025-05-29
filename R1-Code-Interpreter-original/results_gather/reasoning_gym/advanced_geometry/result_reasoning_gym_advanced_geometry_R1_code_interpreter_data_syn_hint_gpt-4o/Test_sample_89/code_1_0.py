# Re-evaluate the equations of the altitudes
eq1 = Eq(y, -7 * (x + 4) + 8)
eq2 = Eq(y, x + 8)

# Solve the equations again
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)