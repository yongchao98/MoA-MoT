# The second voltage plateau of graphite anode occurs at roughly 20% to 50% stoichiometry.
# This corresponds to the phase transition where domains of Stage 2 material (LiC12)
# are formed from a more dilute Lithium-Graphite phase.

# A voltage plateau represents a region of phase coexistence, and its voltage level
# is determined by the chemical potential of the denser phase being formed.
# In this case, the denser phase is Stage 2.

# The chemical potential for intercalation into Stage 2 is denoted by μ_2.

# The relationship between the half-cell voltage (V) and the chemical potential (μ)
# relative to a Li/Li+ reference electrode is given by:
# V = -μ / e
# where 'e' is the elementary charge of a Li+ ion.

# Therefore, substituting μ_2 for the second plateau gives the requested formula.
# The final equation includes the number '2' for the stage number.

stage_number = 2
chemical_potential_symbol = "μ"
charge_symbol = "e"

# Construct and print the formula string for the voltage of the second plateau.
print(f"-{chemical_potential_symbol}_{stage_number} / {charge_symbol}")