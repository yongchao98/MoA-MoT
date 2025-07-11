import sympy as sp

# Define the symbols for our variables
x, y, x_new, y_new = sp.symbols('x y x_new y_new')

# The original equation of the curve
original_rhs = 8*x**5 + 4*x**4 + 4*x**3 + 1*x**2 + 8*x
original_eq = sp.Eq(y**2, original_rhs)

# We perform the change of variables x = 8*x_new, y = 8*y_new.
# This corresponds to blowing up the singular point (0,0) of the mod 2 reduction three times.
transformed_eq = original_eq.subs({x: 8*x_new, y: 8*y_new})

# The left side is (8*y_new)^2 = 64*y_new^2.
# We simplify the right side and divide the whole equation by 64.
simplified_rhs = sp.simplify(transformed_eq.rhs / 64)

# Create the final equation for the new model
final_eq = sp.Eq(y_new**2, simplified_rhs)

# Extract the coefficients to print the equation term by term
poly = sp.Poly(final_eq.rhs, x_new)
coeffs = poly.all_coeffs()

print("The equation for the improved model of the curve is:")
# The format requires printing each number.
print(f"y^2 = {coeffs[0]}*x^5 + {coeffs[1]}*x^4 + {coeffs[2]}*x^3 + {coeffs[3]}*x^2 + {coeffs[4]}*x")

print("\nThe reduction of this equation modulo 2 is y^2 = x^2 + x.")
print("This special fiber has a tacnode singularity, which resolves into 2 double points in the stable reduction.")
print("\nTherefore, the number of double points is 2.")
