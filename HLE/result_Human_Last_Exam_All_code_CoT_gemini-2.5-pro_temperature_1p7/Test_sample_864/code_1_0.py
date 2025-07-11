import numpy as np

# We start by clarifying the assumed model due to issues with the problem statement.
# The given formulas for q_k lead to a non-physical power spectrum since q_1 > q_0.
# We proceed assuming a corrected, well-defined autocorrelation sequence where:
q0 = 5/3
q1 = 2/3
# This gives Q(D) = (2/3)D^{-1} + 5/3 + (2/3)D.

print("--- Step 1: Spectral Factorization of Q(D) ---")
# We want to find a causal and minimum-phase G(D) such that Q(D) = G(D)G(D^{-1}).
# For G(D) = g_0 + g_1*D, we have G(D)G(D^{-1}) = (g_0^2 + g_1^2) + (g_0*g_1)(D + D^{-1}).
# Comparing with Q(D) = 5/3 + (2/3)(D + D^{-1}), we need to solve:
# 1) g_0^2 + g_1^2 = 5/3
# 2) g_0 * g_1 = 2/3
# Solving this system yields g_0^2 = 4/3 and g_1^2 = 1/3 (for the minimum-phase solution).
# Taking the positive roots for simplicity:
g0 = np.sqrt(4/3)
g1 = np.sqrt(1/3)
# The minimum-phase spectral factor is G(D) = sqrt(4/3) + sqrt(1/3)*D.
# This can be simplified to G(D) = (1/sqrt(3)) * (2 + D).

print("The causal, minimum-phase spectral factor found is G(D) = (1/sqrt(3)) * (2 + D).\n")


print("--- Step 2: Determine the Whitening Filter W(D) ---")
# The whitening filter W(D) is chosen to make the output noise white and the
# equivalent channel F(D) = Q(D)W(D) causal.
# The general form is W(D) = D^k / G(D) (with a constant factor set to 1).
# The resulting channel is F(D) = D^k * G(D^{-1}).
# With G(D) being a polynomial of degree L=1, G(D^{-1}) has terms up to D^{-1}.
# To make F(D) causal, we need D^k to eliminate all negative powers of D.
# F(D) = D^k * (g_0 + g_1*D^{-1}) = g_0*D^k + g_1*D^{k-1}.
# The lowest power is k-1, so we need k-1 >= 0, which means k >= 1.
# The simplest choice that ensures causality is k = 1.
k = 1

print(f"The filter form is W(D) = D^k / G(D). To ensure causality, we select k = {k}.")
print("Thus, the whitening filter is W(D) = D / [ (1/sqrt(3)) * (2 + D) ].\n")


print("--- Step 3: Final Whitening Filter Equation ---")
# Simplifying the expression W(D) = D / [ (1/sqrt(3)) * (2 + D) ]
# gives W(D) = (sqrt(3) * D) / (2 + D).
final_g0_coeff = 2.0
final_g1_coeff = 1.0

print("The final expression for the whitening filter is:")
print(f"W(D) = (sqrt(3) * D) / ({final_g0_coeff:.0f} + {final_g1_coeff:.0f} * D)")
print("\nThe numbers in this final equation are:")
print(f"- Numerator coefficient: sqrt(3) ≈ {np.sqrt(3):.4f}")
print(f"- Denominator constant term: {final_g0_coeff:.0f}")
print(f"- Denominator coefficient of D: {final_g1_coeff:.0f}")
