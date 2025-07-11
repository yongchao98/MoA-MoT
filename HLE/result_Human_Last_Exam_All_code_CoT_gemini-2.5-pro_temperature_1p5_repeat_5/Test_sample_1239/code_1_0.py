def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries of the 
    incompressible Navier-Stokes equations in R^3.

    The dimension is the total number of independent generators of the symmetry group.
    """
    
    # Number of generators for each type of symmetry
    time_translation = 1 # t -> t + a
    space_translations = 3 # x -> x + a (vector)
    space_rotations = 3    # x -> Rx (SO(3) rotation matrix)
    galilean_boosts = 3    # x -> x + vt, u -> u + v
    scaling_transformation = 1 # t -> a^2*t, x -> a*x, u -> u/a, p -> p/a^2
    pressure_shift = 1     # p -> p + c
    
    # The total dimension is the sum of the generators from each symmetry subgroup.
    total_dimension = (time_translation + 
                       space_translations + 
                       space_rotations + 
                       galilean_boosts + 
                       scaling_transformation + 
                       pressure_shift)
                       
    print("The dimension of the Lie group is calculated by summing the number of generators for each symmetry:")
    # The final print statement shows each number in the final equation as requested.
    print(f"{time_translation} (time translation) + "
          f"{space_translations} (space translations) + "
          f"{space_rotations} (space rotations) + "
          f"{galilean_boosts} (Galilean boosts) + "
          f"{scaling_transformation} (scaling) + "
          f"{pressure_shift} (pressure shift) = {total_dimension}")

calculate_navier_stokes_symmetry_dimension()
