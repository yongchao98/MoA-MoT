def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    3D incompressible Navier-Stokes equations by summing the dimensions
    of its known symmetry subgroups.
    """
    
    # Dimensions of the various symmetry subgroups
    time_translation = 1
    space_translations = 3
    rotations = 3
    galilean_boosts = 3
    scaling_transformations = 1
    pressure_shifts = 1
    
    symmetries = {
        "Time translations": time_translation,
        "Space translations": space_translations,
        "Rotations": rotations,
        "Galilean boosts": galilean_boosts,
        "Scaling transformations": scaling_transformations,
        "Pressure shifts": pressure_shifts
    }
    
    print("The dimension of the Lie group is the sum of the dimensions of its component symmetries:")
    for name, dim in symmetries.items():
        print(f"- {name}: {dim}")
        
    total_dimension = sum(symmetries.values())
    
    # Create the equation string as requested
    dimension_values = list(symmetries.values())
    equation_string = " + ".join(map(str, dimension_values))
    
    print("\nThe final calculation is:")
    print(f"{equation_string} = {total_dimension}")

if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()
