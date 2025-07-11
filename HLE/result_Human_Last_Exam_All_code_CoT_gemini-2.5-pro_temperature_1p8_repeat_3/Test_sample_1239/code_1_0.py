# This script calculates the dimension of the finite-dimensional Lie group of
# symmetries of the incompressible Navier-Stokes equations in 3D.

# Dimension from time translation (t -> t + a)
dim_time_translation = 1

# Dimension from space translations (x -> x + b) in R^3
dim_space_translation = 3

# Dimension from space rotations (group SO(3))
dim_rotations = 3

# Dimension from Galilean boosts (u -> u + v) in R^3
dim_galilean_boosts = 3

# Dimension from scaling transformations (a single parameter group)
dim_scaling = 1

# The total dimension is the sum of the dimensions of these independent
# symmetry transformations.
total_dimension = (dim_time_translation +
                   dim_space_translation +
                   dim_rotations +
                   dim_galilean_boosts +
                   dim_scaling)

# Print the breakdown and the final sum in an equation format.
print("The total dimension is the sum of dimensions from each symmetry transformation:")
print(f"Time Translation (1) + Space Translations (3) + Rotations (3) + Galilean Boosts (3) + Scaling (1)")
print(f"Equation: {dim_time_translation} + {dim_space_translation} + {dim_rotations} + {dim_galilean_boosts} + {dim_scaling} = {total_dimension}")
