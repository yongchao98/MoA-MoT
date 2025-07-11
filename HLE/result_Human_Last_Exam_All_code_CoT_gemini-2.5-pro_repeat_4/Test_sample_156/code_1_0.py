def find_max_colony_size():
    """
    Calculates the maximum possible area K for a rectangle on an 8x8 board
    with a perimeter of at most 30.

    This is based on the reasoning that:
    1. The colony's perimeter never increases.
    2. The maximum initial perimeter of the 8 starting squares is 30.
    3. Therefore, the final colony's perimeter must be at most 30.
    4. The most area-efficient shape for a given perimeter is a square or rectangle.
    """
    board_dim = 8
    max_perimeter = 30
    max_area = 0
    best_w, best_h = 0, 0

    # Iterate through all possible rectangle dimensions on an 8x8 board
    for w in range(1, board_dim + 1):
        for h in range(1, board_dim + 1):
            perimeter = 2 * w + 2 * h
            
            # Check if the perimeter constraint is met
            if perimeter <= max_perimeter:
                area = w * h
                if area > max_area:
                    max_area = area
                    best_w = w
                    best_h = h
    
    # The maximal size K is the maximum area found.
    # The problem asks to output the numbers in the final equation.
    print(f"The analysis shows the maximal colony size K is bounded by the area of a shape with perimeter <= 30.")
    print(f"The optimal rectangular shape found is a {best_w}x{best_h} rectangle.")
    print(f"Width = {best_w}")
    print(f"Height = {best_h}")
    print(f"Perimeter = 2 * {best_w} + 2 * {best_h} = {2 * best_w + 2 * best_h}")
    print(f"Maximal Area K = {best_w} * {best_h} = {max_area}")

find_max_colony_size()