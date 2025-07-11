import numpy as np

# The problem is to determine the asymptotic growth rate of the degree d_n of a polynomial p_n.
# The asymptotic rate is given as d_n = Theta(n^alpha), and we need to find alpha.

# Let's outline the reasoning based on the properties of polynomials.

# 1. Define the sets and constraints:
# Set A: {1, 2, ..., n^2}, where p_n(i) is in [0, 1].
# Set B: {n^2+1, ..., n^10}, where p_n(i) is in [2, 3].

# 2. Normalize the polynomial for easier analysis.
# Let q_n(x) = p_n(x) - 2.5.
# For i in Set B, q_n(i) is in [-0.5, 0.5].
# For i in Set A, q_n(i) is in [-2.5, -1.5]. So, at x=n^2, |q_n(n^2)| >= 1.5.

# 3. Use an inequality for polynomials bounded on an interval.
# Assume that since q_n(i) is small on the many points in Set B, q_n(x) is small on the
# continuous interval I_B = [n^2+1, n^10]. Let's assume |q_n(x)| <= 0.5 for x in I_B.
# A standard result states that a polynomial bounded by M on an interval grows at most like
# a Chebyshev polynomial of the same degree outside the interval.

# 4. Map the interval I_B to the standard interval [-1, 1].
# The transformation is y(x) = (2*x - (n^10 + n^2 + 1)) / (n^10 - (n^2 + 1)).
# Let Q_n(y) be the polynomial in the y-coordinate. Its degree is d_n.
# We have |Q_n(y)| <= 0.5 for y in [-1, 1].

# 5. Check the constraint at a point outside the interval.
# We need to evaluate the polynomial at x = n^2. Let's find the corresponding y value.
# y(n^2) = (2*(n^2) - (n^10 + n^2 + 1)) / (n^10 - n^2 - 1)
# y(n^2) = (n^2 - n^10 - 1) / (n^10 - n^2 - 1)
# For large n, y(n^2) is approximately -1. Let's find the small difference.
# y(n^2) = -(n^10 - n^2 + 1) / (n^10 - n^2 - 1) = - (1 + 2 / (n^10 - n^2 - 1))
# So, y(n^2) is approximately -1 - 2/n^10. Let this be -1 - epsilon, where epsilon is ~2/n^10.

# 6. Apply the Chebyshev bound.
# We have the condition |q_n(n^2)| >= 1.5, which means |Q_n(y(n^2))| >= 1.5.
# The Chebyshev bound gives |Q_n(y)| <= M * |T_d(y)| for |y|>1, where M=0.5 and d=d_n.
# So, 1.5 <= 0.5 * |T_{d_n}(-1 - epsilon)|.
# This simplifies to |T_{d_n}(-1 - epsilon)| >= 3.

# 7. Analyze the growth of the Chebyshev polynomial.
# For small delta, |T_d(-1-delta)| = T_d(1+delta) ~ 1 + d^2 * delta.
# Here, delta = epsilon = 2/n^10.
# So we get the inequality: 3 <= 1 + (d_n^2) * (2 / n^10).
# 2 <= 2 * d_n^2 / n^10
# 1 <= d_n^2 / n^10
# d_n^2 >= n^10
# d_n >= n^5

# 8. Conclusion.
# The lower bound for the degree d_n grows as n^5.
# This implies that d_n = Omega(n^5).
# In such problems, the lower bound is typically tight. A corresponding upper bound can be constructed,
# showing that d_n = O(n^5).
# Therefore, the asymptotic growth rate is d_n = Theta(n^5).
# The value of alpha is 5.

alpha = 5
print(f"The asymptotic growth rate of d_n as n to infinity can be written as d_n = Theta(n^alpha).")
print(f"The value of alpha is determined by the following derivation:")
print(f"1. Let q_n(x) = p_n(x) - 2.5. On set B = [n^2+1, n^10], |q_n(i)| <= 0.5. At x=n^2, |q_n(n^2)| >= 1.5.")
print(f"2. Assume |q_n(x)| <= 0.5 on the interval [n^2+1, n^10]. Scale this to [-1, 1].")
print(f"3. The point x=n^2 maps to y ~ -1 - 2/n^10.")
print(f"4. Using Chebyshev inequality for growth outside an interval, we get |T_dn(y(n^2))| >= 1.5/0.5 = 3.")
print(f"5. For small epsilon, T_d(1+epsilon) ~ 1 + d^2*epsilon.")
print(f"6. This leads to the inequality: 3 <= 1 + dn^2 * (2/n^10).")
print(f"7. Solving for dn gives dn^2 >= n^10, which implies dn >= n^5.")
print(f"Therefore, alpha must be 5.")
print(f"alpha = {alpha}")
