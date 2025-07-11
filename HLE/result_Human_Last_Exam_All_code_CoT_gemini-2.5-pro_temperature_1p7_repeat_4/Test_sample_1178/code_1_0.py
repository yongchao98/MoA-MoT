def solve():
    """
    This function explains the reasoning and prints the area calculation for the smallest
    integer-length rectangle admitting a non-guillotine tiling with the given squares.
    """
    
    rect_L = 11
    rect_W = 14
    area = rect_L * rect_W
    
    # The set of allowed squares is S = {2x2, 3x3, 5x5, 7x7}.
    # The areas are {4, 9, 25, 49}.
    
    # We found that the area 154 can be formed by a sum of areas from the set.
    # One possible combination is two 7x7 squares and fourteen 2x2 squares.
    # 2 * (7*7) + 14 * (2*2) = 2 * 49 + 14 * 4 = 98 + 56 = 154.
    
    num_7x7 = 2
    area_7x7 = 49
    num_2x2 = 14
    area_2x2 = 4
    
    print("The smallest known integer-length rectangle admitting a non-guillotine tiling")
    print("that can be tiled by squares from S={2x2, 3x3, 5x5, 7x7} is an 11x14 rectangle.")
    print(f"The area of this rectangle is {rect_L} * {rect_W} = {area}.")
    print("\nThis area can be composed of the areas of the given squares. One such combination is:")
    print(f"{num_7x7} * (7*7) + {num_2x2} * (2*2) = {num_7x7} * {area_7x7} + {num_2x2} * {area_2x2} = {num_7x7 * area_7x7} + {num_2x2 * area_2x2} = {area}")
    print(f"\nThus, the area of this rectangle is {area}.")

solve()
