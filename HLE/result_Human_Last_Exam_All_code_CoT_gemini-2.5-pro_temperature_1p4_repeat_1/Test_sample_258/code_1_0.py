def solve_grid_crossing():
    """
    Calculates the minimal and maximal number of grid cells a circle can cross.
    """
    R = 500

    # The number of cells crossed by a circle of integer radius R is given by
    # C = 8R - 4k, where k is the number of half-integers in the set
    # {x0, y0, x0+y0, x0-y0}, where (x0, y0) is the circle's center.

    # To find the maximum number of cells, we need to minimize k.
    # We can choose a center (x0, y0) like (0.1, 0.2) where none of
    # {x0, y0, x0+y0, x0-y0} are half-integers.
    # In this case, k = 0.
    k_min = 0
    max_cells = 8 * R - 4 * k_min

    # To find the minimum number of cells, we need to maximize k.
    # We can choose a center (x0, y0) like (0.5, 0.5).
    # x0 = 0.5 is a half-integer.
    # y0 = 0.5 is a half-integer.
    # x0+y0 = 1 is not a half-integer.
    # x0-y0 = 0 is not a half-integer.
    # The maximum possible value for k is 2.
    k_max = 2
    min_cells = 8 * R - 4 * k_max
    
    # We need to output the equation for both cases.
    print(f"The radius of the circle is R = {R}.")
    print("\nThe number of cells crossed (C) depends on the number of half-integers (k) in {x0, y0, x0+y0, x0-y0}.")
    print("The formula is C = 8*R - 4*k.")
    
    print("\nFor the maximal number of cells, we choose a center where k is minimal (k=0).")
    print(f"Maximal cells = 8 * {R} - 4 * {k_min} = {max_cells}")
    
    print("\nFor the minimal number of cells, we choose a center where k is maximal (k=2).")
    print(f"Minimal cells = 8 * {R} - 4 * {k_max} = {min_cells}")

    print(f"\nThe minimal and maximal numbers of cells are ({min_cells}, {max_cells}).")


solve_grid_crossing()