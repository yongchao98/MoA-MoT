import math

def solve_box_optimization():
    """
    Searches for a more material-efficient box to store energy balls.

    The function iterates through possible integer dimensions for a rectangular box,
    calculates if it can hold at least 27 balls (diameter 4cm), and checks
    if its surface area is less than the original 12x12x12 box (864 cm^2).
    """
    ball_diameter = 4
    min_balls_to_contain = 27
    initial_area = 12 * 12 * 6

    best_box = {
        "l": 12,
        "w": 12,
        "h": 12,
        "area": initial_area,
        "balls": 27
    }

    # Search a reasonable range of dimensions. Very large or elongated
    # boxes will have a much larger surface area.
    # The sum of dimensions for a better box should be similar to or less than the original 12+12+12=36
    search_range = sum(best_box.values()) // 3 * 2 # Heuristic search limit
    search_range = 30 # A reasonable upper bound

    # Using integers for dimensions is a practical starting point for the search.
    # We iterate l <= w <= h to avoid checking permutations of the same box.
    for l in range(1, search_range + 1):
        for w in range(l, search_range + 1):
            # Pruning the search space. For a fixed l,w the area increases with h.
            # If 2*l*w is already > initial_area, no larger w or h can result in a smaller area.
            if 2 * l * w > initial_area:
                break
            for h in range(w, search_range + 1):
                surface_area = 2 * (l*w + w*h + h*l)

                # If the area is already not an improvement, skip to the next
                if surface_area >= best_box["area"]:
                    # Since h is increasing, further iterations in this inner loop won't be better.
                    break
                
                # Calculate how many balls fit in a simple grid packing
                balls_l = math.floor(l / ball_diameter)
                balls_w = math.floor(w / ball_diameter)
                balls_h = math.floor(h / ball_diameter)
                num_balls = balls_l * balls_w * balls_h

                if num_balls >= min_balls_to_contain:
                    # Found a better box
                    best_box["l"] = l
                    best_box["w"] = w
                    best_box["h"] = h
                    best_box["area"] = int(surface_area)
                    best_box["balls"] = num_balls

    if best_box["area"] < initial_area:
        l, w, h, area = best_box["l"], best_box["w"], best_box["h"], best_box["area"]
        print(f"{l}:{w}:{h}:{area}")
    else:
        # After searching, no better box was found.
        # The 12x12x12 cube is optimal for simple grid packing.
        print("0")

solve_box_optimization()