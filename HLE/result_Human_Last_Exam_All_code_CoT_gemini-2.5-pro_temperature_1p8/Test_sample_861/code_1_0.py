import sympy

# Step 1 & 2: Find roots (B, D) and particular solution (E)
# The characteristic equation is 8*r**2 - 6*r + 1 = 0
# The roots are r = (6 +/- sqrt(36-32))/16 = (6 +/- 2)/16, which are 1/2 and 1/4.
# We assign the larger root to B by convention.
B = sympy.Rational(1, 2)
D = sympy.Rational(1, 4)

# The particular solution E is found by solving 8E - 6E + E = 1 => 3E = 1.
E = sympy.Rational(1, 3)

# Step 3: Set up and solve the system of equations for A and C
# The general solution is y[n] = A*B**n + C*D**n + E
# We define A and C as symbolic variables to be solved.
A = sympy.Symbol('A')
C = sympy.Symbol('C')

# Using the initial condition y[0] = 1:
# A * B**0 + C * D**0 + E = 1  => A + C + E = 1
eq1 = sympy.Eq(A + C + E, 1)

# Using the initial condition y[-1] = 2:
# A * B**(-1) + C * D**(-1) + E = 2 => A/B + C/D + E = 2
eq2 = sympy.Eq(A/B + C/D, 2 - E)

# Solve the system for A and C
solution = sympy.solve((eq1, eq2), (A, C))
A_val = solution[A]
C_val = solution[C]

# The solved parameters are now known
# A = A_val, B = B, C = C_val, D = D, E = E

# Print the final form of the difference equation and the parameter values
print("The closed form solution is:")
print(f"y[n] = ({A_val}) * ({B})^n + ({C_val}) * ({D})^n + ({E})")
print("\nWhere the parameters are:")
print(f"A = {A_val}")
print(f"B = {B}")
print(f"C = {C_val}")
print(f"D = {D}")
print(f"E = {E}")

# Step 4: Calculate the final expression
term1 = E / A_val
term2 = (D * C_val) / B
result = term1 + term2

# Print the calculation steps for the final expression
print("\nCalculating the expression E/A + (D*C)/B:")
print(f"E/A = ({E}) / ({A_val}) = {term1}")
print(f"(D*C)/B = (({D})*({C_val})) / ({B}) = {term2}")
print(f"Result = {term1} + {term2} = {result}")
