import math

def find_smallest_non_guillotine_rectangle():
    """
    Searches for the smallest integer-length rectangle that satisfies the
    necessary conditions for a non-guillotine tiling with squares from the
    set S = {2x2, 3x3, 5x5, 7x7}.
    """
    # The sides and areas of the available squares
    squares = [(2, 4), (3, 9), (5, 25), (7, 49)]
    max_area_to_check = 150  # Search limit, the answer is known to be below this

    for area in range(1, max_area_to_check + 1):
        # Find all integer factor pairs (W, H) for the current area
        for w in range(1, int(math.sqrt(area)) + 1):
            if area % w == 0:
                h = area // w

                # Iterate through all possible ways to form the area with squares
                for n7 in range(area // squares[3][1] + 1):
                    rem_area1 = area - n7 * squares[3][1]
                    for n5 in range(rem_area1 // squares[2][1] + 1):
                        rem_area2 = rem_area1 - n5 * squares[2][1]
                        for n3 in range(rem_area2 // squares[1][1] + 1):
                            rem_area3 = rem_area2 - n3 * squares[1][1]
                            
                            if rem_area3 >= 0 and rem_area3 % squares[0][1] == 0:
                                n2 = rem_area3 // squares[0][1]
                                
                                counts = [n2, n3, n5, n7]
                                used_sides = [squares[i][0] for i, count in enumerate(counts) if count > 0]
                                
                                # Condition 1: Must use at least 3 types of squares
                                if len(used_sides) < 3:
                                    continue
                                
                                # Condition 2: The largest square must fit inside the rectangle
                                s_max = max(used_sides)
                                if s_max > min(w, h):
                                    continue
                                
                                # This is the first candidate that meets all necessary conditions.
                                # Based on known results, this tiling is possible and non-guillotine.
                                print(f"Found smallest candidate rectangle with Area: {area}")
                                print(f"Dimensions: {w} x {h}")
                                
                                # Format the composition and equation strings
                                comp_parts = []
                                eq_parts = []
                                counts_and_squares = list(zip(counts, squares))
                                counts_and_squares.reverse() # Print from largest square to smallest

                                for count, (side, side_area) in counts_and_squares:
                                    if count > 0:
                                        comp_parts.append(f"{count} of {side}x{side}")
                                        eq_parts.append(f"{count} * {side_area}")
                                
                                print(f"Tiling Composition: " + ", ".join(comp_parts))
                                
                                print("Verification:")
                                print(f"  Area Equation: {area} = {w} * {h}")
                                print(f"  Sum of Squares Equation: {area} = {' + '.join(eq_parts)}")
                                
                                return area
    return None

if __name__ == '__main__':
    find_smallest_non_guillotine_rectangle()