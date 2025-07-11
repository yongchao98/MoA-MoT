import sympy

# Define x and constants C, D as symbolic variables
x = sympy.Symbol('x')
C = sympy.Symbol('C')
D = sympy.Symbol('D')

# Step 1: State the problem
print("The given differential equation is (dy/dx)⁴ + x(dy/dx) - 3y(x) = 0.")
print("The boundary condition is y(-1) = 0.")
print("-" * 20)

# Step 2: Assume a linear solution y(x) = C*x + D
y = C * x + D
y_prime = sympy.diff(y, x)
print(f"Let's test for a linear solution of the form y(x) = {y}.")
print(f"The derivative is dy/dx = {y_prime}.")
print("-" * 20)

# Step 3: Substitute into the ODE
ode_substituted = y_prime**4 + x * y_prime - 3 * y
print("Substituting into the ODE gives:")
print(f"{y_prime}**4 + x*({y_prime}) - 3*({y}) = 0")
# Simplify the expression
equation = sympy.expand(ode_substituted)
print(f"Simplified equation: {equation} = 0")
print("-" * 20)

# Step 4: For the equation to hold for all x, coefficients must be zero
# Extract coefficient of x
coeff_x = equation.coeff(x)
print(f"The coefficient of x is {coeff_x}.")
print(f"Setting it to zero: {coeff_x} = 0, which implies C = 0.")
C_val = 0

# Substitute C=0 to find D
constant_term = equation.subs(x, 0)
print(f"The constant term is {constant_term}.")
constant_term_eq = constant_term.subs(C, C_val)
print(f"With C = {C_val}, we get: {constant_term_eq} = 0, which implies D = 0.")
D_val = 0
print("-" * 20)

# Step 5: State the solution
print(f"The only linear solution is y(x) = {C_val}*x + {D_val}, which is y(x) = 0.")
final_y_at_0 = 0
print(f"This solution satisfies the boundary condition y(-1) = 0.")
print("-" * 20)

# Step 6: Final answer
print("Therefore, the membrane's deflection at x = 0 is:")
print(f"y(0) = {final_y_at_0}")