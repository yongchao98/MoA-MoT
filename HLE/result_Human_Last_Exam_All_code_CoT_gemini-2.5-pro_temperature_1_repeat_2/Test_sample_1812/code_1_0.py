# The user wants to compute the determinant of the given matrix.
# The matrix is:
# A = [[ 2, -1, -1],
#      [-1,  2, -1],
#      [-1, -1,  2]]
#
# The determinant of a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]] is
# calculated as a(ei - fh) - b(di - fg) + c(dh - eg).
# I will implement this formula and print out the steps.

# Define the matrix elements
a, b, c = 2, -1, -1
d, e, f = -1, 2, -1
g, h, i = -1, -1, 2

# Calculate the determinant
determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

# Print the calculation steps
print("To compute the determinant of the matrix:")
print(f"A = [[{a}, {b}, {c}],")
print(f"     [{d}, {e}, {f}],")
print(f"     [{g}, {h}, {i}]]")
print("")

print("The formula for the determinant is: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
print("")

print("Substituting the values from matrix A:")
print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))")
print("")

print("Calculating the terms inside the parentheses:")
print(f"det(A) = {a}({e*i} - {f*h}) - ({b})({d*i} - {f*g}) + ({c})({d*h} - {e*g})")
print(f"det(A) = {a}({e*i - f*h}) - ({b})({d*i - f*g}) + ({c})({d*h - e*g})")
print("")

print("Performing the multiplication:")
print(f"det(A) = {a * (e*i - f*h)} - ({-b * (d*i - f*g)}) + ({c * (d*h - e*g)})")
print("")

print("Finally, summing the terms:")
print(f"det(A) = {a * (e*i - f*h) - (b * (d*i - f*g)) + (c * (d*h - e*g))}")
print("")

print(f"The determinant of the matrix is: {determinant}")
<<<0>>>