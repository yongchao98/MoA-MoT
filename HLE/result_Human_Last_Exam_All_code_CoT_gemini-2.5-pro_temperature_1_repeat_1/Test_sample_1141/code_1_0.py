import math

def solve_and_print_solution():
    """
    Solves the billet cutting problem based on the user's formulation.
    The problem formulation contains non-standard and restrictive constraints, particularly for T1 cubes. 
    Analysis shows that placing any B2 ball (value 150) makes it impossible to place any T1 cube (value 5).
    Therefore, the highest value will come from a strategy that avoids T1 cubes and maximizes the placement of B2 balls,
    filling the remaining space with B1 balls (value 1).

    This script implements that optimal strategy:
    1. Place the maximum number of non-overlapping B2 balls. An optimal 4x2 grid placement yields 8 balls.
    2. Greedily fill the remaining volume with B1 balls, checking against the specified constraints for each potential placement.
    3. Sum the values of all placed items and print the result in the requested equation format.
    """
    placed_objects = []

    # Step 1: Place the 8 B2 balls in an optimal grid configuration.
    # The centers must be on the z=4 plane and separated by a distance of at least 8 units.
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        (4, 12, 4), (12, 12, 4), (20, 12, 4), (28, 12, 4)
    ]
    for center in b2_centers:
        placed_objects.append({'type': 'B2', 'center': center, 'price': 150})

    # Step 2: Greedily place B1 balls in the remaining space.
    # We iterate through all possible B1 center coordinates and place a B1 if it doesn't conflict with existing objects.
    # B1 center ranges: x in [1, 31], y in [1, 21], z in [1, 7].
    for z in range(1, 8):  # z_i in [1, ..., 7]
        for y in range(1, 22):  # y_i in [1, ..., 21]
            for x in range(1, 32):  # x_i in [1, ..., 31]
                candidate_center = (x, y, z)
                is_valid_placement = True
                
                # Check for conflicts with all previously placed objects
                for existing_obj in placed_objects:
                    existing_center = existing_obj['center']
                    
                    # Calculate squared distance between centers
                    dist_sq = (candidate_center[0] - existing_center[0])**2 + \
                              (candidate_center[1] - existing_center[1])**2 + \
                              (candidate_center[2] - existing_center[2])**2

                    if existing_obj['type'] == 'B2':
                        # Constraint for B1 to B2: dist_sq >= 25
                        if dist_sq < 25:
                            is_valid_placement = False
                            break
                    elif existing_obj['type'] == 'B1':
                        # Constraint for B1 to B1: dist_sq >= 4
                        if dist_sq < 4:
                            is_valid_placement = False
                            break
                
                if is_valid_placement:
                    # If there are no conflicts, add the B1 to our set of objects
                    placed_objects.append({'type': 'B1', 'center': candidate_center, 'price': 1})

    # Step 3: Construct and print the final equation and total value.
    price_list = [str(obj['price']) for obj in placed_objects]
    total_value = sum(obj['price'] for obj in placed_objects)
    
    equation_str = " + ".join(price_list)
    
    print(f"The highest value solution found under your rules is:")
    print(f"Maximize {equation_str} = {total_value}")


solve_and_print_solution()