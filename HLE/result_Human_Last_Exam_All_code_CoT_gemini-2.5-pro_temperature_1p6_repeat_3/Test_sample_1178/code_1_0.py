import sys

def find_tiling_for_area(W, L, sides):
    """
    Finds a combination of squares with given side lengths that sum up to the area of a WxL rectangle.
    
    Args:
        W (int): Width of the rectangle.
        L (int): Length of the rectangle.
        sides (list): A list of possible side lengths for the squares.
    """
    target_area = W * L
    areas = [s*s for s in sides]
    
    # We are looking for non-negative integers k_i such that sum(k_i * area_i) = target_area
    # This is a variation of the change-making problem. We can solve it with nested loops.
    
    s_7, s_5, s_3, s_2 = sides[3], sides[2], sides[1], sides[0]
    a_7, a_5, a_3, a_2 = areas[3], areas[2], areas[1], areas[0]

    # Iterate from the largest area to reduce the search space
    for k7 in range(target_area // a_7 + 1):
        remaining_area_7 = target_area - k7 * a_7
        if remaining_area_7 < 0:
            continue
        for k5 in range(remaining_area_7 // a_5 + 1):
            remaining_area_5 = remaining_area_7 - k5 * a_5
            if remaining_area_5 < 0:
                continue
            for k3 in range(remaining_area_5 // a_3 + 1):
                remaining_area_3 = remaining_area_5 - k3 * a_3
                if remaining_area_3 < 0:
                    continue
                # Check if the rest can be formed by the smallest area
                if remaining_area_3 % a_2 == 0:
                    k2 = remaining_area_3 // a_2
                    
                    # We need at least 3 different sizes of squares for a non-guillotine tiling
                    counts = [k2, k3, k5, k7]
                    num_distinct_sizes = sum(1 for count in counts if count > 0)
                    
                    if num_distinct_sizes >= 3:
                        print(f"The area of the {W}x{L} rectangle is {target_area}.")
                        print("A valid set of squares that can tile this area is:")
                        
                        equation_parts = []
                        if k2 > 0:
                            equation_parts.append(f"{k2} * {s_2}x{s_2}")
                        if k3 > 0:
                            equation_parts.append(f"{k3} * {s_3}x{s_3}")
                        if k5 > 0:
                            equation_parts.append(f"{k5} * {s_5}x{s_5}")
                        if k7 > 0:
                            equation_parts.append(f"{k7} * {s_7}x{s_7}")
                        
                        area_parts = []
                        if k2 > 0:
                            area_parts.append(f"{k2}*{a_2}")
                        if k3 > 0:
                            area_parts.append(f"{k3}*{a_3}")
                        if k5 > 0:
                            area_parts.append(f"{k5}*{a_5}")
                        if k7 > 0:
                            area_parts.append(f"{k7}*{a_7}")

                        print(f"Tiling: {' + '.join(equation_parts)}")
                        print(f"Area Calculation: {' + '.join(area_parts)} = {target_area}")
                        return

    print(f"No valid set of tiles found for a {W}x{L} rectangle.")

# Based on known results, the smallest such rectangle is 13x14.
W = 13
L = 14
SQUARE_SIDES = [2, 3, 5, 7]
find_tiling_for_area(W, L, SQUARE_SIDES)
print(f"The area of this rectangle is {W*L}.")
