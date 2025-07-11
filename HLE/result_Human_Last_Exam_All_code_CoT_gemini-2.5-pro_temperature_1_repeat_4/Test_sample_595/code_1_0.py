def solve_triangle_grid_problem():
    """
    This script calculates the largest number k of coordinate grid squares
    that the perimeter of a triangle with side lengths 18, 18, and 18*sqrt(2)
    can pass through without containing any lattice points.

    The solution is based on placing the triangle's legs parallel to the
    coordinate axes to maximize grid crossings and then using the
    Principle of Inclusion-Exclusion to count the unique squares.
    """

    # The triangle is a right isosceles triangle. We place it with legs
    # parallel to the axes, shifted by a small epsilon to avoid lattice points.
    # Let vertices be B(0.1, 0.1), A(18.1, 0.1), C(0.1, 18.1).

    # Number of squares crossed by the horizontal leg AB (length 18)
    # It passes through squares with y-index 0, and x-indices from 0 to 18.
    num_squares_AB = 18 + 1
    
    # Number of squares crossed by the vertical leg BC (length 18)
    # It passes through squares with x-index 0, and y-indices from 0 to 18.
    num_squares_BC = 18 + 1

    # Number of squares crossed by the hypotenuse AC
    # It crosses 18 vertical grid lines and 18 horizontal grid lines.
    num_squares_AC = 18 + 18

    # Overlap between AB and BC is the single square containing vertex B: (0,0).
    overlap_AB_BC = 1

    # Overlap between AB and AC. The hypotenuse starts in square (18,0) and
    # passes through (17,0) before crossing the y=1 line. Both are on AB's path.
    overlap_AB_AC = 2

    # Overlap between BC and AC. Symmetrically, the hypotenuse ends in square (0,18)
    # after passing through (0,17). Both are on BC's path.
    overlap_BC_AC = 2

    # The three-way overlap is zero, as the hypotenuse doesn't pass through (0,0).
    overlap_ABC = 0

    # Calculate the total number of unique squares using the Principle of Inclusion-Exclusion.
    total_squares = (num_squares_AB + num_squares_BC + num_squares_AC -
                     (overlap_AB_BC + overlap_AB_AC + overlap_BC_AC) +
                     overlap_ABC)

    print("Calculation of the largest number of squares, k:")
    print("-------------------------------------------------")
    print(f"Squares crossed by leg AB: {num_squares_AB}")
    print(f"Squares crossed by leg BC: {num_squares_BC}")
    print(f"Squares crossed by hypotenuse AC: {num_squares_AC}")
    print(f"Overlap of squares for AB and BC: {overlap_AB_BC}")
    print(f"Overlap of squares for AB and AC: {overlap_AB_AC}")
    print(f"Overlap of squares for BC and AC: {overlap_BC_AC}")
    print(f"Overlap of squares for AB, BC, and AC: {overlap_ABC}")
    print("\nApplying the Principle of Inclusion-Exclusion:")
    print(f"k = |AB| + |BC| + |AC| - (|AB∩BC| + |AB∩AC| + |BC∩AC|) + |AB∩BC∩AC|")
    print(f"k = {num_squares_AB} + {num_squares_BC} + {num_squares_AC} - ({overlap_AB_BC} + {overlap_AB_AC} + {overlap_BC_AC}) + {overlap_ABC}")
    print(f"k = {num_squares_AB + num_squares_BC + num_squares_AC} - {overlap_AB_BC + overlap_AB_AC + overlap_BC_AC} + {overlap_ABC}")
    print(f"k = {total_squares}")

solve_triangle_grid_problem()
<<<69>>>