import sys

def solve_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the 3D incompressible Navier-Stokes equations.
    """
    print("This program determines the dimension of the Lie group of symmetries of the 3D incompressible Navier-Stokes equations.")
    print("The dimension is the sum of the number of independent generators for each type of symmetry transformation.")
    print("-" * 50)

    # Dimension from time translation
    dim_time = 1
    print(f"Symmetry: Time translation (t -> t + c)")
    print(f"Dimension: {dim_time}\n")

    # Dimension from space translations
    dim_space = 3
    print(f"Symmetry: Space translations (x -> x + c) in R^3")
    print(f"Dimension: {dim_space}\n")

    # Dimension from space rotations
    dim_rotation = 3
    print(f"Symmetry: Space rotations (SO(3) group)")
    print(f"Dimension: {dim_rotation}\n")

    # Dimension from Galilean boosts
    dim_galilean = 3
    print(f"Symmetry: Galilean boosts (change of inertial frame)")
    print(f"Dimension: {dim_galilean}\n")

    # Dimension from scaling
    dim_scaling = 1
    print(f"Symmetry: Scaling or dilation")
    print(f"Dimension: {dim_scaling}\n")
    
    # Dimension from pressure shift
    dim_pressure = 1
    print(f"Symmetry: Constant pressure shift (p -> p + c)")
    print(f"Reason: The equations only depend on the pressure gradient ∇p.")
    print(f"Dimension: {dim_pressure}\n")

    # Summing up all dimensions
    total_dimension = dim_time + dim_space + dim_rotation + dim_galilean + dim_scaling + dim_pressure
    
    print("-" * 50)
    print("The total dimension is the sum of the dimensions from each independent symmetry.")
    
    # The final print statement that shows the equation explicitly as requested.
    print(f"Total Dimension = {dim_time} + {dim_space} + {dim_rotation} + {dim_galilean} + {dim_scaling} + {dim_pressure} = {total_dimension}")

if __name__ == '__main__':
    solve_navier_stokes_symmetry_dimension()