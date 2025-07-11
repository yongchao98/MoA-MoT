import sys

def solve_lie_group_dimension():
    """
    Calculates the dimension of the principal Lie group of symmetries
    of the incompressible Navier-Stokes equations in R^3.
    """

    # Introduction explaining the context
    print("The symmetry group of the incompressible Navier-Stokes equations in R^3 is, strictly speaking, infinite-dimensional due to pressure gauge freedom.", file=sys.stdout)
    print("However, the dimension of its maximal finite-dimensional subalgebra (the principal Lie algebra) is a well-defined integer value.", file=sys.stdout)
    print("\nThis dimension is the sum of the number of generators from different types of symmetries:", file=sys.stdout)

    # Number of generators for each type of symmetry
    time_translation = 1
    space_translations = 3
    spatial_rotations = 3
    galilean_boosts = 3
    scaling = 1
    pressure_shift = 1 # Constant pressure shift

    # Calculate the total dimension
    total_dimension = (time_translation + space_translations + spatial_rotations +
                       galilean_boosts + scaling + pressure_shift)

    # Print the detailed breakdown of the calculation
    print(f"\nDimension = (Time Translations) + (Space Translations) + (Spatial Rotations) + (Galilean Boosts) + (Scaling) + (Constant Pressure Shift)", file=sys.stdout)
    print(f"Dimension = {time_translation} + {space_translations} + {spatial_rotations} + {galilean_boosts} + {scaling} + {pressure_shift}", file=sys.stdout)
    print(f"Total Dimension = {total_dimension}", file=sys.stdout)

solve_lie_group_dimension()