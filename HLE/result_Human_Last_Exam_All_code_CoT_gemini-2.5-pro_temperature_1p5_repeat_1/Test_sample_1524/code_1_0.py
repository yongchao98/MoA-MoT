# The average number of bosons in a state with energy ε is n(ε) = 1 / (exp((ε - μ) / (k_B*T)) - 1).
# For n(ε) to be a positive real number, the denominator must be positive:
# exp((ε - μ) / (k_B*T)) - 1 > 0
# This implies that the argument of the exponential must be positive:
# ε - μ > 0  or  μ < ε

# This must hold for all energy levels ε. The most restrictive case is the lowest energy level,
# the ground state energy, which we denote as ε_0.

# Let's assign a value to the ground state energy for this example.
# We can set it to 0, or a small positive value. Let's use a symbolic name.
ground_state_energy_symbol = "ε_0"
chemical_potential_symbol = "μ"

# The final derived relationship is that the chemical potential must be less than the ground state energy.
# At the limit (during condensation or at T=0), it can be equal.
print("The fundamental limit for the chemical potential (μ) is derived from the Bose-Einstein distribution.")
print("The occupation number of any state must be non-negative, which leads to the inequality:")
print(f"    {chemical_potential_symbol} < ε")
print("\nThis must be true for ALL energy states ε. The most restrictive condition is set by the lowest possible energy, the ground state energy ε_0.")
print("Therefore, the fundamental limit is:")
print(f"    {chemical_potential_symbol} <= {ground_state_energy_symbol}")

# At absolute zero (T=0), all particles in a non-interacting Bose gas occupy the ground state.
# At this temperature, the chemical potential is exactly equal to the ground state energy.
# μ(T=0) = ε_0

print("\nAt zero temperature, the chemical potential is exactly equal to the ground state energy: μ(T=0) = ε_0.")
print("Thus, the fundamental limit on the chemical potential is that it must be less than or equal to the chemical potential of a non-interacting Bose gas at zero temperature.")
