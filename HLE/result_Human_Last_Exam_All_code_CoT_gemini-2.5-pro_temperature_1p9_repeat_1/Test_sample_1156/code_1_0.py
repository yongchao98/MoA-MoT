import sympy

# Define symbols
x, C = sympy.symbols('x C')

# Define the unnormalized density function form
# From literature, the form is C * (1/sqrt(x) - 1/2)
rho_unnormalized = C * (x**(-1/2) - 1/2)

# To find C, we need to solve the normalization integral equation:
# Integral from 0 to 1 of rho(x) dx = 1
integral_value = sympy.integrate(rho_unnormalized / C, (x, 0, 1))

# The equation is C * integral_value = 1. We solve for C.
C_val = 1 / integral_value

# Get the numerator and denominator for the coefficients in the final equation
a = (C_val * rho_unnormalized.coeff(x**(-1/2))).simplify()
b = (C_val * rho_unnormalized.coeff(x, 0)).simplify() # coeff for the constant term

# Output the equation for the density
# To display fractions nicely, we get numer/denom from SymPy objects
a_num, a_den = sympy.fraction(a)
b_num, b_den = sympy.fraction(b)

print("The normalised density of the invariant measure is rho(x) = a/sqrt(x) + b")
print(f"where the coefficients are a = {a} and b = {b}.")
print("\nThe final equation is:")
print(f"rho(x) = ({a_num}/{a_den})/sqrt(x) + ({b_num}/{b_den})")
print(f"rho(x) = {a_num}/({a_den}*sqrt(x)) - {abs(b_num)}/{b_den}")