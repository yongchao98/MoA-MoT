# The dimension of the Lie group of symmetries of the incompressible
# Navier-Stokes equations in R^3 is the sum of the dimensions of its
# constituent symmetry subgroups.

# Dimension of the time translation group (t' = t + a)
time_translations_dim = 1

# Dimension of the space translation group (x' = x + a) in R^3
space_translations_dim = 3

# Dimension of the spatial rotation group SO(3)
spatial_rotations_dim = 3

# Dimension of the Galilean boosts group (v' = v + v_0) in R^3
galilean_boosts_dim = 3

# Dimension of the scaling transformations group
scaling_transformations_dim = 1

# Dimension of the pressure shift group (p' = p + c)
pressure_shift_dim = 1

# The total dimension is the sum of these individual dimensions.
total_dimension = (time_translations_dim + 
                   space_translations_dim + 
                   spatial_rotations_dim + 
                   galilean_boosts_dim + 
                   scaling_transformations_dim + 
                   pressure_shift_dim)

# Display the calculation, showing each number in the final equation.
print("The total dimension is the sum of the dimensions from each symmetry group:")
print(f"Total Dimension = (Time Translations) + (Space Translations) + (Spatial Rotations) + (Galilean Boosts) + (Scaling) + (Pressure Shift)")
print(f"Total Dimension = {time_translations_dim} + {space_translations_dim} + {spatial_rotations_dim} + {galilean_boosts_dim} + {scaling_transformations_dim} + {pressure_shift_dim}")
print(f"\nThe final calculated dimension is: {total_dimension}")

print("\n<<<12>>>")