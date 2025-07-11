# The dimension of the Lie group of symmetries of the incompressible Navier-Stokes equations
# in R^3 is the sum of the dimensions of its independent one-parameter subgroups.

# Number of parameters for time translation (t -> t + a)
time_translations = 1

# Number of parameters for space translations (x -> x + b) in R^3
space_translations = 3

# Number of parameters for spatial rotations in R^3 (group SO(3))
spatial_rotations = 3

# Number of parameters for Galilean boosts (x -> x - v*t) in R^3
galilean_boosts = 3

# Number of parameters for the pressure shift (p -> p + c)
pressure_shift = 1

# Number of parameters for the scaling transformation of space, time, velocity, and pressure
scaling_symmetry = 1

# The total dimension is the sum of the dimensions of all these symmetries.
total_dimension = (time_translations + 
                   space_translations + 
                   spatial_rotations + 
                   galilean_boosts + 
                   pressure_shift + 
                   scaling_symmetry)

# Print the full equation as requested
print("The dimension of the Lie group of symmetries is the sum of the dimensions of each independent transformation:")
print(f"{time_translations} (time translation) + "
      f"{space_translations} (space translations) + "
      f"{spatial_rotations} (spatial rotations) + "
      f"{galilean_boosts} (Galilean boosts) + "
      f"{pressure_shift} (pressure shift) + "
      f"{scaling_symmetry} (scaling)")
print(f"= {total_dimension}")