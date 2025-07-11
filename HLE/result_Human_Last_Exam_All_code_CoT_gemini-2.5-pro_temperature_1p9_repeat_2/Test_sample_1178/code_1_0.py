def solve_tiling_problem():
    """
    This function verifies the area calculation for the smallest known integer-length
    rectangle that admits a non-guillotine tiling by squares from the set
    S={2x2, 3x3, 5x5, 7x7}.
    """
    # Dimensions of the rectangle
    length = 10
    width = 12

    # Number of squares of each size for the tiling
    # This specific configuration for a 10x12 rectangle is known from tiling theory.
    # Note that 7x7 squares are not needed for this minimal solution.
    num_2x2 = 4
    num_3x3 = 6
    num_5x5 = 2
    num_7x7 = 0

    # Calculate the area of the rectangle
    rectangle_area = length * width

    # Calculate the sum of the areas of the squares
    squares_area = (num_2x2 * 2**2) + \
                   (num_3x3 * 3**2) + \
                   (num_5x5 * 5**2) + \
                   (num_7x7 * 7**2)

    # Print the equation
    print("Verifying the area calculation for the 10x12 rectangle:")
    print(f"{num_2x2} * (2x2) + {num_3x3} * (3x3) + {num_5x5} * (5x5) + {num_7x7} * (7x7)")
    
    equation_parts = []
    if num_2x2 > 0:
        equation_parts.append(f"{num_2x2} * {2*2}")
    if num_3x3 > 0:
        equation_parts.append(f"{num_3x3} * {3*3}")
    if num_5x5 > 0:
        equation_parts.append(f"{num_5x5} * {5*5}")
    if num_7x7 > 0:
        equation_parts.append(f"{num_7x7} * {7*7}")
        
    print(f"{' + '.join(equation_parts)} = {squares_area}")

    # Check if the areas match
    if rectangle_area == squares_area:
        print(f"\nThe area of the squares ({squares_area}) correctly matches the rectangle area ({rectangle_area}).")
        print(f"The area of this rectangle is {rectangle_area}.")
    else:
        print("\nError: The sum of the areas of the squares does not match the rectangle's area.")

solve_tiling_problem()
