import math
from fractions import Fraction

# Step 1: Find the particular solution E
# The equation is 8y[n] - 6y[n-1] + y[n-2] = 1.
# Assume a particular solution y_p[n] = E (a constant).
# 8*E - 6*E + E = 1
# 3*E = 1 => E = 1/3
E = Fraction(1, 3)
print(f"Step 1: Find the particular solution.")
print(f"Assuming a constant solution y_p[n] = E, we get 8E - 6E + E = 1, so 3E = 1.")
print(f"The value of E is {E}\n")

# Step 2: Find the homogeneous solution roots B and D
# The characteristic equation is 8r^2 - 6r + 1 = 0.
# We solve for r using the quadratic formula: r = (-b ± sqrt(b^2 - 4ac)) / 2a
a, b, c = 8, -6, 1
discriminant = b**2 - 4*a*c
r1 = (-b + math.sqrt(discriminant)) / (2*a)
r2 = (-b - math.sqrt(discriminant)) / (2*a)

# By convention, assign the larger root to B
B = Fraction(r1) if r1 > r2 else Fraction(r2)
D = Fraction(r2) if r1 > r2 else Fraction(r1)
print(f"Step 2: Find the roots of the characteristic equation 8r^2 - 6r + 1 = 0.")
print(f"The roots are r1 = {B} and r2 = {D}.")
print(f"Assigning the larger root to B, we have B = {B} and D = {D}.\n")

# Step 3: Use initial conditions to find A and C
# The general solution is y[n] = A*(B)^n + C*(D)^n + E
# Initial conditions: y[0] = 1 and y[-1] = 2.
# For n=0: y[0] = A*(B)^0 + C*(D)^0 + E = 1
# A + C + E = 1 => A + C = 1 - E
# For n=-1: y[-1] = A*(B)^-1 + C*(D)^-1 + E = 2
# A/B + C/D + E = 2 => A/B + C/D = 2 - E

# We have a system of two linear equations for A and C:
# 1) A + C = 1 - E
# 2) A/B + C/D = 2 - E
eq1_rhs = 1 - E
eq2_rhs = 2 - E

# From eq1: A = eq1_rhs - C
# Substitute into eq2: (eq1_rhs - C)/B + C/D = eq2_rhs
# eq1_rhs/B - C/B + C/D = eq2_rhs
# C*(1/D - 1/B) = eq2_rhs - eq1_rhs/B
# C = (eq2_rhs - eq1_rhs/B) / (1/D - 1/B)

C = (eq2_rhs - eq1_rhs / B) / (1/D - 1/B)
A = eq1_rhs - C

print(f"Step 3: Solve for A and C using initial conditions y[0]=1 and y[-1]=2.")
print(f"Equation from y[0]=1: A + C = 1 - E = {eq1_rhs}")
print(f"Equation from y[-1]=2: A/B + C/D = 2 - E = {eq2_rhs}")
print(f"Solving this system gives A = {A} and C = {C}.\n")

# The final equation is:
print("The final closed form solution is:")
print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")
print(f"where A={A}, B={B}, C={C}, D={D}, E={E}\n")


# Step 4: Calculate the final expression
# Expression to calculate: E/A + (D*C)/B
result = E/A + (D*C)/B
print("Step 4: Calculate the final expression E/A + (D*C)/B.")
print(f"E/A = ({E}) / ({A}) = {E/A}")
print(f"(D*C)/B = (({D}) * ({C})) / ({B}) = {(D*C)/B}")
print(f"The result of E/A + (D*C)/B is {result}")
print(f"The result as a decimal is {float(result)}")
