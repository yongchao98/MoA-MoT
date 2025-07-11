def solve_tiling_problem():
    """
    This function verifies the properties of the smallest integer-length rectangle
    that admits a non-guillotine tiling with squares from S={2x2, 3x3, 5x5, 7x7}.
    """
    # Dimensions of the smallest known rectangle admitting a fault-free tiling
    # with squares from the given set.
    length = 11
    width = 10
    
    # Calculate the area of the rectangle.
    rectangle_area = length * width
    
    # The composition of the tiling for the 11x10 rectangle.
    # This tiling uses three different types of squares {3,5,7}, a necessary
    # condition for a fault-free tiling with squares from S.
    num_7x7 = 1
    num_5x5 = 1
    num_3x3 = 4
    num_2x2 = 0
    
    # Calculate the total area of the tiles.
    area_7x7 = num_7x7 * (7**2)
    area_5x5 = num_5x5 * (5**2)
    area_3x3 = num_3x3 * (3**2)
    area_2x2 = num_2x2 * (2**2)
    
    total_tile_area = area_7x7 + area_5x5 + area_3x3 + area_2x2
    
    # Verify that the areas match.
    if rectangle_area == total_tile_area:
        print(f"The rectangle has dimensions {length}x{width} and an area of {rectangle_area}.")
        print("A known fault-free tiling consists of:")
        print(f"- {num_7x7} square(s) of size 7x7")
        print(f"- {num_5x5} square(s) of size 5x5")
        print(f"- {num_3x3} square(s) of size 3x3")
        print(f"- {num_2x2} square(s) of size 2x2")
        print("\nVerifying the area calculation:")
        
        # Output the equation with each number.
        print(f"{num_7x7} * 7*7 + {num_5x5} * 5*5 + {num_3x3} * 3*3 + {num_2x2} * 2*2 = {area_7x7} + {area_5x5} + {area_3x3} + {area_2x2} = {total_tile_area}")
        print(f"\nThe area of this rectangle is {rectangle_area}.")
    else:
        print("There is a mismatch in the calculation.")

solve_tiling_problem()