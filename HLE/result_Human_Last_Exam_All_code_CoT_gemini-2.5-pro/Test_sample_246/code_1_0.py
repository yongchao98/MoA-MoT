# The manifold is the double branched cover of S^4 over the 5-twist-spun trefoil.
# We want the minimal number of generators of its fundamental group, d(pi_1(M)).
# This is bounded by the number of generators of the first homology group, H_1(M).
# The order of H_1(M) is given by |Delta_K(-1)|, where Delta_K is the Alexander polynomial of the knot.

# The knot is the 5-twist-spun trefoil, T_5.
# Its Alexander polynomial is Delta_T5(t) = Delta_T(t) * (t^5 - 1) / (t - 1).
# The Alexander polynomial of the trefoil is Delta_T(t) = t - 1 + 1/t.

# We evaluate at t = -1.
t = -1
n = 5

# Calculate Delta_T(-1)
delta_T_at_minus_1 = t - 1 + 1/t
print(f"The Alexander polynomial of the trefoil, Delta_T(t), evaluated at t = -1 is: {int(delta_T_at_minus_1)}")

# Calculate the factor (t^n - 1) / (t - 1) at t = -1
factor_at_minus_1 = (t**n - 1) / (t - 1)
print(f"The factor (t^{n}-1)/(t-1) for n={n} evaluated at t = -1 is: {int(factor_at_minus_1)}")

# Calculate Delta_T5(-1)
delta_T5_at_minus_1 = delta_T_at_minus_1 * factor_at_minus_1
print(f"The Alexander polynomial of the 5-twist-spun trefoil, Delta_T5(t), evaluated at t = -1 is: {int(delta_T5_at_minus_1)}")

# The order of the first homology group H_1(M) is the absolute value.
order_H1 = abs(delta_T5_at_minus_1)
print(f"The order of the first homology group H_1(M) is: {int(order_H1)}")

# Since |H_1(M)| = 3, H_1(M) must be the cyclic group Z_3.
# The fundamental group pi_1(M) is a quotient of pi_1(Sigma_2(T)) = Z_3.
# Since H_1 is non-trivial, pi_1(M) must also be Z_3.
# The minimal number of generators for the cyclic group Z_3 is 1.
min_generators = 1
print(f"The minimal number of generators of the fundamental group is: {min_generators}")