# The open-circuit voltage (V) of the half-cell is related to the chemical potential
# of the reference electrode (μ_ref) and the electrode itself (μ_electrode) by:
# e * V = μ_ref - μ_electrode
#
# A voltage plateau corresponds to a phase transition where μ_electrode is constant.
#
# The problem defines μ_k as the chemical potential of intercalating at stage k.
#
# The first plateau (V1 ≈ 0.09V) is the transition from Stage 2 to Stage 1.
# This involves intercalating Li+ into a Stage 2 structure.
# So, for this plateau, μ_electrode = μ_2.
# e * V1 = μ_ref - μ_2  (1)
#
# The second plateau (V2 ≈ 0.13V) is the transition from Stage 3 to Stage 2.
# This involves intercalating Li+ into a Stage 3 structure.
# So, for this plateau, μ_electrode = μ_3.
# e * V2 = μ_ref - μ_3  (2)
#
# To find a formula for V2 that involves different μ_k, we can eliminate μ_ref.
# From (1), we have: μ_ref = e * V1 + μ_2
# Substitute this into (2):
# e * V2 = (e * V1 + μ_2) - μ_3
#
# Divide by e:
# V2 = V1 + (μ_2 - μ_3) / e
#
# Substitute the given value for V1 (0.09V):
# V2 = 0.09 + (μ_2 - μ_3) / e
#
# This formula expresses the voltage of the second plateau in terms of the voltage
# of the first plateau, the chemical potentials for intercalating at stages 2 and 3,
# and the elementary charge.

# The final formula will be printed as a string.
final_formula = "0.09 + (μ_2 - μ_3) / e"

print(final_formula)