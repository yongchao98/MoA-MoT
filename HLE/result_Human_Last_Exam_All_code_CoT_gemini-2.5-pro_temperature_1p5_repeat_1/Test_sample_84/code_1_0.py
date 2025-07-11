import math

# The problem is to determine the asymptotic growth rate alpha for the minimal degree d_n
# of a polynomial p_n satisfying certain conditions on two sets of points.
# The sets are S1 = {1, 2, ..., n^2} and S2 = {n^2+1, ..., n^10}.
# The conditions are p_n(i) in [0, 1] for i in S1, and p_n(i) in [2, 3] for i in S2.

# Plan:
# 1. Establish an upper bound for d_n by constructing a suitable polynomial.
#    This construction is based on Chebyshev polynomials and shows d_n = O(n).
#    This implies alpha <= 1.
# 2. Establish a lower bound for d_n.
#    This is done by showing that a polynomial of degree o(n) cannot bridge the gap
#    between the two sets of values fast enough.
#    This proves d_n = Omega(n), which implies alpha >= 1.
# 3. Combining the upper and lower bounds gives d_n = Theta(n).
# 4. From d_n = Theta(n^alpha), we conclude that alpha = 1.

# Upper bound derivation:
# Let's construct a polynomial p_n(x) using Chebyshev polynomial T_d(y).
# Map the interval [1, n^2] to [-1, 1] using y(x) = (2x - n^2 - 1) / (n^2 - 1).
# Let p_n(x) = 0.5 * T_d(y(x)) + 0.5.
# For i in S1, y(i) is in [-1, 1], so p_n(i) is in [0, 1].
# For i in S2, we need p_n(i) in [2, 3], which means T_d(y(i)) in [3, 5].
# The tightest constraint is at i = n^2 + 1.
# y(n^2 + 1) = (n^2 + 1) / (n^2 - 1) = 1 + 2 / (n^2 - 1).
# For large n, T_d(1 + epsilon) approx 0.5 * exp(d * sqrt(2*epsilon)).
# Here, epsilon is approx 2/n^2.
# We need 0.5 * T_d(y(n^2+1)) + 0.5 >= 2 => T_d(y(n^2+1)) >= 3.
# 0.5 * exp(d * sqrt(4/n^2)) >= 3 => exp(2*d/n) >= 6 => d >= n * ln(6)/2.
# So, a degree of d_n = O(n) is sufficient. Thus, alpha <= 1.

# Lower bound derivation:
# Assume d_n = o(n).
# Let p_n be a polynomial of degree d_n.
# p_n is bounded by 1 on the n^2 integer points in [1, n^2].
# Because d_n << n^2, p_n is bounded on the whole interval [1, n^2] by a constant K.
# By polynomial growth estimates, for x > n^2:
# |p_n(x)| <= K * |T_{d_n}(y(x))|.
# At x = n^2 + 1, we require p_n(n^2+1) >= 2.
# 2 <= K * T_{d_n}(y(n^2+1)) = K * T_{d_n}(1 + 2/(n^2-1)).
# 2 <= K * 0.5 * exp(d_n * sqrt(4/n^2)) = (K/2) * exp(2*d_n/n).
# ln(4/K) <= 2*d_n/n => d_n >= (n/2) * ln(4/K).
# This means d_n = Omega(n), which contradicts d_n = o(n).
# Therefore, alpha >= 1.

# Conclusion:
# Since alpha <= 1 and alpha >= 1, we must have alpha = 1.

alpha = 1.0

print(f"The asymptotic growth rate is d_n = Theta(n^alpha).")
print(f"The value of alpha is {alpha}.")
