import fractions

# The coefficients of the polynomial f(x) = x^3 + ax^2 + bx + c were derived
# from the problem's conditions.
a = fractions.Fraction(-3, 8)
b = fractions.Fraction(-5, 8)
c = fractions.Fraction(0, 1)

# We need to compute f(3).
x = 3
f_3 = x**3 + a * x**2 + b * x + c

# The final equation is f(3) = 3^3 + (-3/8)*3^2 + (-5/8)*3 + 0
val_x3 = x**3
val_x2 = x**2
val_x = x

# Print the final calculation step-by-step as requested
print(f"f(3) = {val_x3} + ({a})*({val_x2}) + ({b})*({val_x}) + {c} = {f_3}")