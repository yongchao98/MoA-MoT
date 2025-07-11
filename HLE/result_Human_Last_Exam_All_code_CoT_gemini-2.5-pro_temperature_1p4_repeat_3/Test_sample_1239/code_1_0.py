#
# Calculate the dimension of the Lie group of symmetries for the
# 3D incompressible Navier-Stokes equations.
#
# The total dimension is the sum of the dimensions of the individual
# symmetry subgroups.
#

# Dimension from time translations (t' = t + a)
dim_time_translation = 1

# Dimension from space translations (x' = x + a) in R^3
dim_space_translation = 3

# Dimension from spatial rotations (x' = Rx) in R^3
dim_rotation = 3

# Dimension from Galilean boosts (x' = x + vt) in R^3
dim_galilean_boost = 3

# Dimension from scaling transformations
dim_scaling = 1

# Calculate the total dimension of the Lie group
total_dimension = (
    dim_time_translation +
    dim_space_translation +
    dim_rotation +
    dim_galilean_boost +
    dim_scaling
)

# Print the breakdown of the calculation and the final result.
# The final equation shows each number being added.
print("The dimension is the sum of dimensions from:")
print("- Time Translations:   ", dim_time_translation)
print("- Space Translations:  ", dim_space_translation)
print("- Spatial Rotations:   ", dim_rotation)
print("- Galilean Boosts:     ", dim_galilean_boost)
print("- Scaling:             ", dim_scaling)
print("\nFinal Equation:")
print(f"{dim_time_translation} + {dim_space_translation} + {dim_rotation} + {dim_galilean_boost} + {dim_scaling} = {total_dimension}")