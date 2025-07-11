import math

# The problem reduces to finding the largest alpha that satisfies either E_2(alpha) = 0 or psi_2(alpha, alpha) = 0.
# A quasi-exactly-solvable model analysis leads to two conditions for alpha.

# Condition 1: E_2(alpha) = 0. This leads to the equation 5*alpha^2 - 8 = 0.
# The numbers in this equation are 5, 8, 0.
a1 = 5
b1 = 0
c1 = -8
# Solve the quadratic equation a1*alpha^2 + c1 = 0
# alpha^2 = -c1 / a1
alpha_sq_1 = -c1 / a1
alpha1 = math.sqrt(alpha_sq_1)

# Condition 2: psi_2(alpha, alpha) = 0. This leads to the equation 2*alpha^4 - 2*alpha^3 - 1 = 0.
# We can find the root of this polynomial numerically.
# Let's define the polynomial function.
def poly(a):
    return 2 * a**4 - 2 * a**3 - 1

# Simple numerical root finding (bisection method)
low = 1.0
high = 2.0
while high - low > 1e-9:
    mid = (low + high) / 2
    if poly(mid) < 0:
        low = mid
    else:
        high = mid
alpha2 = (low + high) / 2

# The problem asks for the largest value, alpha_0.
alpha_0 = max(alpha1, alpha2)

# The question asks to output each number in the final equation.
# The largest alpha comes from the first equation: 5*alpha^2 - 8 = 0.
# Let's print the numbers from that equation and the final answer.

print(f"The equation for the largest alpha_0 is: {a1}*alpha^2 - {abs(c1)} = 0")
print(f"The largest value alpha_0 such that F(alpha_0) = 0 is: {alpha_0}")
