def calculate_origami_regions():
    """
    Calculates the number of regions the paper is divided into by the
    creases of a standard origami crane.
    """

    # Step 1: Define the number of lines (creases) for a standard origami crane.
    # The standard crane crease pattern has:
    # - 2 main diagonal lines
    # - 2 median lines (horizontal and vertical)
    # - 4 lines forming a central rotated square (a "diamond")
    # - 8 short lines connecting the vertices of the central diamond to the
    #   edges of the paper.
    # Total lines (L) = 2 + 2 + 4 + 8 = 16
    num_lines = 16

    # Step 2: Define the number of intersection points for these creases.
    # The intersection points are:
    # - 1 central point where the 2 diagonals and 2 medians meet.
    # - 4 points where the diagonals intersect the sides of the central diamond.
    # - 4 points where the medians intersect the sides of the central diamond.
    # Total points (P) = 1 + 4 + 4 = 9
    num_points = 9

    # Step 3: Calculate the number of regions using Euler's formula for planar graphs,
    # which simplifies to R = 1 + L + P for this type of problem.
    num_regions = 1 + num_lines + num_points

    # Step 4: Print the breakdown and the final calculation.
    print("To find the number of regions, we use the formula: R = 1 + L + P")
    print(f"Where L is the number of lines and P is the number of intersection points.")
    print("-" * 30)
    print(f"Number of crease lines (L) on the unfolded paper: {num_lines}")
    print(f"Number of intersection points (P) on the unfolded paper: {num_points}")
    print("-" * 30)
    print(f"The equation is: Number of Regions = 1 + {num_lines} + {num_points}")
    print(f"Total number of regions = {num_regions}")

calculate_origami_regions()