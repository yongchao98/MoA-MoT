import sympy as sp

# Define the variable of integration
x = sp.Symbol('x')

# Define the initial condition u(x,y,z,0)
u_initial = - (3 * (2 * sp.exp(x) + 1) * sp.exp(x + sp.Symbol('y') + sp.Symbol('z'))) / \
            ((sp.exp(x) + 1) * sp.exp(x + sp.Symbol('y') + sp.Symbol('z')) + 1)

# Substitute z = -y into the initial condition to get the integrand
# This step simplifies the expression significantly
integrand = u_initial.subs(sp.Symbol('z'), -sp.Symbol('y'))

print("The problem is to compute the integral of u(x,y,-y,1) from x=0 to x=1.")
print("Based on the structure of the initial condition, we assume the solution is stationary on the slice z=-y.")
print("Thus, u(x,y,-y,1) = u(x,y,-y,0).")
print("\nThe integral to compute is:")
sp.pprint(sp.Integral(integrand, (x, 0, 1)))
print("\n")

# Perform the integration symbolically
result = sp.integrate(integrand, (x, 0, 1))

# The result can be expressed in terms of logarithms.
# Let's get the components of the result A*log(B) + C*log(D)
# Sympy might give a complicated form, let's simplify it.
# The simplified form is 3*ln(3) - 3*ln(e^2 + e + 1)
# or 3*ln(3 / (e^2 + e + 1))
term1_coeff = 3
term1_arg = 3
term2_coeff = -3
term2_arg = sp.exp(2) + sp.exp(1) + 1

print("The exact result of the integration is:")
# We manually construct the simplified pretty print for clarity
final_expr = term1_coeff * sp.log(term1_arg) + term2_coeff * sp.log(sp.exp(2)+sp.exp(1)+1)
sp.pprint(final_expr)
print("\nThis can also be written as:")
sp.pprint(term1_coeff * sp.log(term1_arg / term2_arg))


# To output each number as requested:
print("\nThe numbers in the final equation 3*ln(3) - 3*ln(exp(2)+exp(1)+1) are:")
print(f"Coefficient 1: {abs(term1_coeff)}")
print(f"Argument 1: {term1_arg}")
print(f"Coefficient 2: {abs(term2_coeff)}")
print("Argument 2 involves powers of e: exp(2)+exp(1)+1")


# For the final answer, we provide the numerical value.
numerical_result = result.evalf()
print(f"\nThe numerical value of the integral is approximately: {numerical_result}")
