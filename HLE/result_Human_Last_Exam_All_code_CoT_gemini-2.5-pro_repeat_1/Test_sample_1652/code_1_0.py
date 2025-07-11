# This script calculates the initial speed 'u' using Newton's method.
# The problem is modeled by the quadratic equation A*u^2 + B*u + C = 0.
#
# Based on the physics:
# A = sin(120)
# B = 2 * v * sin(60)
# C = -d * g
#
# Using rational approximations:
# g = 9.8 -> 49/5
# sin(60) = sin(120) -> 13/15
# And given values v=5, d=300.
#
# The equation becomes:
# (13/15)u^2 + (2*5*13/15)u - (300*49/5) = 0
# (13/15)u^2 + (26/3)u - 2940 = 0
#
# Multiplying by 15 to clear denominators gives the simplified integer equation:
# 13*u^2 + 130*u - 44100 = 0

# The coefficients of the simplified equation A*u^2 + B*u + C = 0 are:
A = 13
B = 130
C = -44100

print("The simplified equation to solve is A*u^2 + B*u + C = 0, with coefficients:")
print(f"A: {A}")
print(f"B: {B}")
print(f"C: {C}")

# We solve for u using Newton's method: u_new = u_old - f(u) / f'(u)
# where f(u) = 13*u^2 + 130*u - 44100
# and f'(u) = 26*u + 130

# Initial guess for the speed u
u = 50.0

# Iterate 5 times for convergence
for _ in range(5):
    f_u = A * u**2 + B * u + C
    f_prime_u = 2 * A * u + B
    # Avoid division by zero, though not expected here
    if f_prime_u == 0:
        break
    u = u - f_u / f_prime_u

# The problem asks for the result rounded to one decimal place.
u_rounded = round(u, 1)

print("\nResult:")
print(f"The initial speed u is approximately: {u_rounded} m/s")