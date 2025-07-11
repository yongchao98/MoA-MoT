def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    incompressible Navier-Stokes equations in a given number of spatial dimensions.
    """
    # The problem is specified for 3 spatial dimensions.
    n = 3

    # 1. Time translations: t' = t + s
    # This transformation has one parameter (s).
    dim_time_translations = 1

    # 2. Space translations: x' = x + a
    # This transformation has n parameters (the components of vector a).
    dim_space_translations = n

    # 3. Space rotations: x' = R*x, u' = R*u
    # The group of rotations in n dimensions, SO(n), has dimension n*(n-1)/2.
    dim_space_rotations = n * (n - 1) // 2

    # 4. Galilean boosts: x' = x + v*t, u' = u + v
    # This transformation has n parameters (the components of vector v).
    dim_galilean_boosts = n

    # 5. Scaling transformations: t' = a^2*t, x' = a*x, u' = a^-1*u, p' = a^-2*p
    # This transformation has one parameter (a).
    dim_scaling = 1
    
    # 6. Constant pressure shifts: p' = p + c
    # The equations are invariant if a constant is added to the pressure,
    # as only the gradient of pressure appears. This has one parameter (c).
    dim_pressure_shift = 1

    # The total dimension is the sum of the dimensions of all independent symmetries.
    total_dimension = (dim_time_translations +
                       dim_space_translations +
                       dim_space_rotations +
                       dim_galilean_boosts +
                       dim_scaling +
                       dim_pressure_shift)

    print(f"The dimension of the Lie group of symmetries for the incompressible Navier-Stokes equations in {n}D is calculated as follows:")
    print("-" * 50)
    print(f"Dimension from time translations: {dim_time_translations}")
    print(f"Dimension from space translations: {dim_space_translations}")
    print(f"Dimension from space rotations: {dim_space_rotations}")
    print(f"Dimension from Galilean boosts: {dim_galilean_boosts}")
    print(f"Dimension from scaling transformations: {dim_scaling}")
    print(f"Dimension from constant pressure shifts: {dim_pressure_shift}")
    print("-" * 50)
    print("The final equation for the total dimension is:")
    print(f"{dim_time_translations} + {dim_space_translations} + {dim_space_rotations} + {dim_galilean_boosts} + {dim_scaling} + {dim_pressure_shift} = {total_dimension}")

if __name__ == '__main__':
    calculate_navier_stokes_symmetry_dimension()