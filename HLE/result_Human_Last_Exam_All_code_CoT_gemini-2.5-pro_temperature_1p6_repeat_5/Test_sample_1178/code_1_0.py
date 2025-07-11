def find_area_decomposition():
    """
    This function searches for a combination of squares from the set
    S = {2x2, 3x3, 5x5, 7x7} whose areas sum up to the target area of 169,
    which corresponds to a 13x13 rectangle.
    """
    target_area = 169
    squares = [
        {'side': 2, 'area': 4},
        {'side': 3, 'area': 9},
        {'side': 5, 'area': 25},
        {'side': 7, 'area': 49}
    ]

    # We iterate from the largest squares to find a solution
    area_7 = squares[3]['area']
    area_5 = squares[2]['area']
    area_3 = squares[1]['area']
    area_2 = squares[0]['area']

    for n7 in range(target_area // area_7, -1, -1):
        remaining_1 = target_area - n7 * area_7
        for n5 in range(remaining_1 // area_5, -1, -1):
            remaining_2 = remaining_1 - n5 * area_5
            for n3 in range(remaining_2 // area_3, -1, -1):
                remaining_3 = remaining_2 - n3 * area_3
                if remaining_3 % area_2 == 0:
                    n2 = remaining_3 // area_2
                    
                    print(f"The area of the smallest rectangle is {target_area}.")
                    print("A possible decomposition of this area into the given square areas is:")
                    
                    # Store the terms of the equation
                    terms = []
                    if n7 > 0:
                        terms.append(f"{n7} * {area_7}")
                    if n5 > 0:
                        terms.append(f"{n5} * {area_5}")
                    if n3 > 0:
                        terms.append(f"{n3} * {area_3}")
                    if n2 > 0:
                        terms.append(f"{n2} * {area_2}")

                    # Print the equation
                    print(f"{target_area} = {' + '.join(terms)}")
                    return

find_area_decomposition()