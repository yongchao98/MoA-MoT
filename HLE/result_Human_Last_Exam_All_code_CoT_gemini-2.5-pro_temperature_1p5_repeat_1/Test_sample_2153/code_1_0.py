# Values determined from the cheminformatics analysis as described in the plan.

# The molecule identified is Cytidine.
# The Hosoya Z index (Z) for Cytidine is the total number of matchings in its molecular graph.
hosoya_z_index = 353

# The Zagreb(1) index (M1) for Cytidine is the sum of the squares of the degrees of its non-hydrogen atoms.
zagreb1_index = 88

# The problem asks for twice the Hosoya Z to the Zagreb(1) index ratio.
# Formula: 2 * Z / M1
ratio = (2 * hosoya_z_index) / zagreb1_index

# Final output in the requested format
print(f"2 * Hosoya Z / Zagreb(1) for Cytidine:")
print(f"2 * {hosoya_z_index} / {zagreb1_index} = {ratio}")