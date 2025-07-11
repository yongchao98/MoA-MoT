# 1. Define the number of NV orientations based on the lattice type.
# In a hypothetical cubic lattice, the NV axes align with the three
# orthogonal crystal axes <100>.
num_orientations = 3

# 2. Apply an E-field along one axis (e.g., z-axis). We then categorize
# the orientations relative to this field.

# The number of NV orientations PERPENDICULAR to the E-field.
# If E is along [001], the axes [100] and [010] are perpendicular.
# A perpendicular E-field splits the resonance into a doublet (2 lines).
# Both orientations are equivalent and produce the same doublet.
num_resonances_from_perpendicular_group = 2

# The number of NV orientations PARALLEL to the E-field.
# If E is along [001], the axis [001] is parallel.
# A parallel E-field does not split the resonance, only shifts it,
# resulting in a singlet (1 line).
num_resonances_from_parallel_group = 1

# 3. The total number of unique resonances is the sum of lines from
# the distinct groups (the doublet from the perpendicular set and the
# singlet from the parallel set).
total_resonances = num_resonances_from_perpendicular_group + num_resonances_from_parallel_group

# 4. Print the final calculation and result.
print("Calculation for the total number of ODMR resonances:")
print(f"Number of lines from the NV group perpendicular to E-field: {num_resonances_from_perpendicular_group}")
print(f"Number of lines from the NV group parallel to E-field: {num_resonances_from_parallel_group}")

print("\nThe final equation is:")
print(f"{num_resonances_from_perpendicular_group} + {num_resonances_from_parallel_group} = {total_resonances}")