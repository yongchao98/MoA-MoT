# This script will derive and print the formula for the second voltage plateau.

# The formula involves the chemical potentials for forming stage 2 and stage 3.
# Let's define the numbers representing these stages.
stage_k1 = 2
stage_k2 = 3

# The formula is an average, so we will be dividing by 2.
divisor = 2

# We construct the formula as a string.
# V_II: Voltage of the second plateau.
# μ_Li: Chemical potential of the Li/Li+ reference electrode.
# μ_k: Chemical potential for intercalating at stage k.
# e: Elementary charge of a Li+ ion.
# The formula represents the voltage as determined by the reference potential
# minus the average of the chemical potentials for the two relevant transitions,
# all divided by the elementary charge.
final_formula = f"V_II = (μ_Li - (μ_{stage_k1} + μ_{stage_k2}) / {divisor}) / e"

# Print the final formula.
print(final_formula)