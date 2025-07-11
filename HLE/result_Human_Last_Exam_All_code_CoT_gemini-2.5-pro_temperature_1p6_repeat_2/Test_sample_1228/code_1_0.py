# This script calculates the neutrino-antineutrino asymmetry from the decay of long-lived neutral kaons (K_L).
# This asymmetry is a direct consequence of CP violation in the neutral kaon system.

# The asymmetry arises because the long-lived kaon, K_L, is an unequal mixture of K⁰ and K̅⁰.
# The semileptonic decay rate of K_L to neutrinos is slightly different from its decay rate to antineutrinos.

# The experimentally measured value for the real part of the CP violation parameter in K-mixing is denoted Re(ε).
Re_epsilon = 1.66e-3

# The semileptonic charge asymmetry, δ_L, is defined as:
# δ_L = [Γ(K_L → π⁻l⁺ν) - Γ(K_L → π⁺l⁻ν̅)] / [Γ(K_L → π⁻l⁺ν) + Γ(K_L → π⁺l⁻ν̅)]
# where Γ represents the decay rate. This is related to Re(ε) by the formula: δ_L ≈ 2 * Re(ε).

# We will now calculate the value of δ_L.
delta_L = 2 * Re_epsilon

print("The asymmetry between neutrinos and antineutrinos from K_L decays is calculated below.")
print("This asymmetry is quantified by the semileptonic charge asymmetry, δ_L.")
print("\nFormula: δ_L ≈ 2 * Re(ε)")
print("Where Re(ε) is the real part of the CP-violation parameter in neutral kaon mixing.")
print("\nCalculation:")
print(f"δ_L ≈ {2} * {Re_epsilon}")
print(f"δ_L ≈ {delta_L:.5f}")
print(f"\nThis result shows a positive asymmetry of approximately {delta_L*100:.3f}%, indicating a small excess of neutrinos produced from K_L decays.")
