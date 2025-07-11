import sympy

# Define symbols for the equation
rho_prime = sympy.Symbol("rho/rho_0")
x_prime = sympy.Symbol("x/L")
M0 = sympy.sqrt(3)
gamma = sympy.Rational(5, 3)

# Calculate the dimensionless coefficient for the argument of tanh
# Pr = 3/4
# L/delta = (gamma/Pr) / ( (2*(4/3)*(gamma+1)) / (3*(1/2)**2) ) is not quite right
# From derivation: coeff = (L/delta_shock)
# L = (mu*gamma)/(Pr*rho0*u0)
# delta_shock = (2 * (4/3)*mu * (gamma+1)*u0) / ((rho0+rho1)*(u0-u1)**2)
# rho1/rho0 = 2, u1/u0 = 1/2
# rho0+rho1 = 3*rho0
# u0-u1 = u0/2
# delta_shock = (2 * (4/3)*mu * (8/3)*u0) / (3*rho0 * (u0/2)**2) = (256*mu)/(27*rho0*u0)
# L/delta_shock = (mu*(5/3)/((3/4)*rho0*u0)) / ((256*mu)/(27*rho0*u0))
# L/delta_shock = ((20/9)*mu/(rho0*u0)) / ((256/27)*mu/(rho0*u0))
# L/delta_shock = (20/9) * (27/256) = (20*3)/256 = 60/256 = 15/64
coeff = sympy.Rational(15, 64)

# Velocity profile u/u0 = u_prime = 1 / rho_prime
# u_prime = u1/u0 + (1-u1/u0)/(1+exp(2*x/delta)) = 1/2 + (1/2)/(1+exp(2*(15/64)*x_prime))
# ... which simplifies to the tanh form

# Construct the equation from the simplified density profile derived in the plan.
# rho_prime = 4 / (3 - tanh( (15/64) * x_prime ))
equation = sympy.Eq(rho_prime, 4 / (3 - sympy.tanh(coeff * x_prime)))

# Print the final analytical solution
# We explicitly construct the string to ensure the numbers are presented as derived.
# Using sympy.pretty just formats the object, but we want to show the specific numbers.
# The `repr` function can get a string representation from sympy.
lhs = "rho/rho_0"
rhs_numerator = "4"
rhs_denominator = f"3 - tanh(({coeff.p}/{coeff.q}) * x/L)"
final_equation_str = f"{lhs} = {rhs_numerator} / ({rhs_denominator})"

print("The analytical solution for the density profile is:")
print(final_equation_str)

# To be extra clear and fulfill the "output each number" requirement, let's print them one by one.
print("\nFinal equation with explicit numbers:")
print(f"rho/rho_0 = {4} / ({3} - tanh(({15}/{64}) * x/L))")