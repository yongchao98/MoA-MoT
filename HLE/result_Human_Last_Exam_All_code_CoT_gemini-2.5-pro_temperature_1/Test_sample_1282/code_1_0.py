import math

# This script calculates the finite blow-up time for a simplified
# version of the given partial differential equation.
# The analysis of a simplified linear model, \partial_t v = -(1+t) \Delta v,
# shows that a solution with a smooth initial condition v(x, 0) = exp(-a * x^2)
# will blow up in finite time.
#
# The blow-up time T is the smallest positive root of the quadratic equation:
# T^2 + 2*T - 1/(2*a) = 0.
#
# We will solve this equation for a specific parameter choice, a = 0.5.

# Parameter 'a' from the initial condition v(x,0) = exp(-a * x^2)
a = 0.5

# Coefficients of the quadratic equation T^2 + 2T - 1/(2a) = 0
# The equation is c1*T^2 + c2*T + c3 = 0
c1 = 1.0
c2 = 2.0
c3 = -1.0 / (2.0 * a)

print("The blow-up time T is found by solving the following quadratic equation for T > 0:")
print(f"({c1}) * T^2 + ({c2}) * T + ({c3}) = 0")

# We use the quadratic formula to find the roots:
# T = (-c2 +/- sqrt(c2^2 - 4*c1*c3)) / (2*c1)
discriminant = c2**2 - 4 * c1 * c3

if discriminant < 0:
    print("\nNo real solution for T exists.")
else:
    # We are looking for a positive time t > 0
    t1 = (-c2 + math.sqrt(discriminant)) / (2 * c1)
    t2 = (-c2 - math.sqrt(discriminant)) / (2 * c1)

    positive_roots = []
    if t1 > 0:
        positive_roots.append(t1)
    if t2 > 0:
        positive_roots.append(t2)

    if not positive_roots:
        print("\nNo positive finite blow-up time found.")
    else:
        blowup_time = min(positive_roots)
        print(f"\nThe finite blow-up time is T = {blowup_time}")