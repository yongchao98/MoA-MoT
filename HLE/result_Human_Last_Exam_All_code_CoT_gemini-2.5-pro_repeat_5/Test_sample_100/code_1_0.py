import sympy

# Define pi as a symbolic variable
pi = sympy.pi

# Step 1: Define the results of the four definite integrals derived from the main problem.
# The main integral I is split into four parts: I = I1 + I2 + I3 + I4

# I1 corresponds to the term p / (exp(p) - 1)
# The definite integral from 0 to infinity is Gamma(2)*zeta(2) = pi^2/6
I1 = pi**2 / 6

# I2 corresponds to the term p^7 / (exp(p) - 1)
# The definite integral is Gamma(8)*zeta(8) = 7! * (pi^8/9450) = 8*pi^8/15
I2 = (sympy.S(8)/15) * pi**8

# I3 corresponds to the term p*exp(-p) / (exp(p) - 1)
# This integral evaluates to pi^2/6 - 1
I3 = pi**2 / 6 - 1

# I4 corresponds to the term (exp(p/4) - exp(-p/4)) / (2*(exp(p) - 1))
# This integral evaluates to 2 - pi/2
I4 = 2 - pi / 2

# Step 2: Sum the results of the four integrals to get the final answer.
total_value = I1 + I2 + I3 + I4

# Step 3: Simplify the symbolic expression.
simplified_value = sympy.simplify(total_value)

# Step 4: Output the components of the final equation as requested.
# The simplified expression is of the form c1*pi^8 + c2*pi^2 + c3*pi + c4
c_pi8 = simplified_value.coeff(pi**8)
c_pi2 = simplified_value.coeff(pi**2)
c_pi1 = simplified_value.coeff(pi)
constant_term = simplified_value.as_coeff_add(pi)[0]


print("The final result is an expression of the form: c1*pi^8 + c2*pi^2 + c3*pi + c4")
print("The components of the final equation are:")
print(f"c1 (coefficient of pi^8) = {c_pi8}")
print(f"c2 (coefficient of pi^2) = {c_pi2}")
print(f"c3 (coefficient of pi) = {c_pi1}")
print(f"c4 (constant term) = {constant_term}")
print("\nSo the final equation is:")
print(f"{c_pi8}*pi**8 + {c_pi2}*pi**2 + ({c_pi1})*pi + {constant_term}")
