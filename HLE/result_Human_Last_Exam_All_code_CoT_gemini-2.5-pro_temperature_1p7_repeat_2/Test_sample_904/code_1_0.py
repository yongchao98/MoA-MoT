import sympy as sp

# Define symbols
gamma, r = sp.symbols('gamma r')
xi = sp.Function('xi')(r)
Delta_P = sp.Symbol('Delta_P')

# Define the derivatives
xi_prime = xi.diff(r)
xi_double_prime = xi.diff(r, r)

# The target form of the equation is A(r)*xi'' + B(r)*xi' + C(r, xi) = 0
# From our derivation, the equation is:
# gamma * xi'' + (gamma/r) * xi' - Delta_P = 0

# Identify the coefficients by comparison
A_r = gamma
B_r = gamma / r
C_r_xi = -Delta_P

# The final equation is constructed from these parts
final_equation = A_r * xi_double_prime + B_r * xi_prime + C_r_xi

print("The governing linear equation for the interfacial shape xi(r) is derived from the linearized Young-Laplace equation.")
print("The equation has the form: A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")

print("The derived equation is:")
# The sp.pretty_print function provides a more readable output for the equation
sp.pretty_print(sp.Eq(final_equation, 0))

print("\nBy comparing the derived equation to the target form, we can identify the coefficients A(r) and B(r).\n")

print(f"The coefficient of the second derivative, A(r), is: {A_r}")
print(f"The coefficient of the first derivative, B(r), is: {B_r}")
