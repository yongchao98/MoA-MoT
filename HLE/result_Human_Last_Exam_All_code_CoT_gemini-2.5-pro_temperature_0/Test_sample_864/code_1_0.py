import math

# Coefficients for the numerator and denominator of the filter W(D)
# W(D) = c / (1 - a * D^{-1})
c = math.sqrt(3) / 2
a = 1/2

# Print the equation for the whitening filter W(D)
# The equation is derived from the spectral factorization of the assumed Q(D).
# W(D) = F_max(D)^{-1} which makes the equivalent channel H_eq(D) = F_min(D) causal.
print("Based on the corrected channel autocorrelation Q(D) = 5/3 - 2/3*(D + D^{-1}), the whitening filter W(D) is:")
print(f"W(D) = (sqrt(3)/2) / (1 - (1/2)*D^{-1})")
print(f"     = {c:.4f} / (1 - {a:.1f}*D^{-1})")
