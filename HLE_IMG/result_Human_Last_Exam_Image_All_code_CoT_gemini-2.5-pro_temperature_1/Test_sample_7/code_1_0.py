# The second voltage plateau in a graphite anode occurs at roughly
# 20-50% state of charge, which corresponds to the electrochemical
# transition where the material changes from Stage 3 (LiC18) to Stage 2 (LiC12).

# The voltage (V) of the anode half-cell is determined by the
# chemical potential (μ) of Lithium in the graphite, according to the
# formula V = -μ/e, where e is the elementary charge of a Li+ ion.

# A voltage plateau indicates a region where two phases coexist in equilibrium.
# For the second plateau, these phases are Stage 2 and Stage 3. A simple
# and common approximation for the chemical potential during this transition
# is the arithmetic average of the potentials of the two pure stages.

# Let μ_k be the chemical potential of Li+ in Stage k.
# The chemical potential for the transition between Stage 2 and Stage 3 is:
# μ_plateau ≈ (μ_2 + μ_3) / 2

# Therefore, the voltage V for the second plateau is approximated by the formula below.

# Define the integer numbers that appear in the final formula.
stage_number_1 = 2
stage_number_2 = 3
denominator = 2

# Print the final formula clearly, showing each number.
print("The simple formula that best approximates the second plateau voltage (V) is:")
print(f"V = - (μ_{stage_number_1} + μ_{stage_number_2}) / ({denominator} * e)")
