import numpy as np

# Based on the corrected problem where Q(D) = 10/3 + D + D^-1
# We found G(D) = g0 + g1*D
g0_sq = 3.0
g1_sq = 1.0/3.0
g0 = np.sqrt(g0_sq)
g1 = np.sqrt(g1_sq)

# G(D) = g0 + g1*D
# W(D) = 1 / G(D^-1) = 1 / (g0 + g1*D^-1) = D / (g0*D + g1)
# W(D) = (1/g0) * D / (D + g1/g0)
# Numerator of W(D) is [1/g0, 0] corresponding to (1/g0)*D
# Denominator of W(D) is [1, g1/g0] corresponding to D + g1/g0

# Let's write W(D) = (c1*D + c0) / (d1*D + d0)
# From W(D) = D / (g0*D + g1)
c1 = 1
c0 = 0
d1 = g0
d0 = g1

print("Based on the corrected problem where Q(D) = 10/3 + D + D^-1,")
print("the whitening filter W(D) is found to be a rational function in D.")
print("W(D) = (c1*D + c0) / (d1*D + d0)")
print(f"c1 = {c1}")
print(f"c0 = {c0}")
print(f"d1 = {d1:.4f}")
print(f"d0 = {d0:.4f}")

# Final expression: W(D) = D / (sqrt(3)*D + 1/sqrt(3))
# Let's write it out clearly.
print("\nThe whitening filter is W(D) = N(D)/M(D) where:")
print(f"N(D) = {c1:.4f} * D + {c0:.4f}")
print(f"M(D) = {d1:.4f} * D + {d0:.4f}")
print(f"So, W(D) = D / ({d1:.4f}*D + {d0:.4f})")