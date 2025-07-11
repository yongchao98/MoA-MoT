import math

def solve_cutting_problem():
    """
    Calculates the maximum value by placing B2 balls first, then filling the
    rest of the space with B1 balls, according to the specified constraints.
    """
    
    # Step 1: Define the optimal placement for the 12 B2 balls.
    # A pre-calculated dense packing configuration for B2 balls (radius 4)
    # in the allowed volume. z-center is fixed at 4.
    b2_centers = []
    for x in [4, 12, 20, 28]:
        for y in [4, 11, 18]:
            b2_centers.append({'x': x, 'y': y, 'z': 4, 'type': 'B2', 'price': 150})
    
    num_b2 = len(b2_centers)
    b2_value = num_b2 * 150

    # Step 2: Greedily place B1 balls in the remaining space.
    b1_centers = []
    # Iterate through all possible B1 center locations in a fixed order.
    # B1 radius = 1 unit.
    for z in range(1, 8):
        for y in range(1, 22):
            for x in range(1, 32):
                p_b1 = {'x': x, 'y': y, 'z': z}
                
                # Check for collision with B2 balls.
                # Constraint: Squared distance >= (radius_b1 + radius_b2)^2 = (1+4)^2 = 25
                is_valid = True
                for p_b2 in b2_centers:
                    dist_sq = (p_b1['x'] - p_b2['x'])**2 + (p_b1['y'] - p_b2['y'])**2 + (p_b1['z'] - p_b2['z'])**2
                    if dist_sq < 25:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # Check for collision with already placed B1 balls.
                # Constraint: Squared distance >= (radius_b1 + radius_b1)^2 = (1+1)^2 = 4
                for p_b1_placed in b1_centers:
                    dist_sq = (p_b1['x'] - p_b1_placed['x'])**2 + (p_b1['y'] - p_b1_placed['y'])**2 + (p_b1['z'] - p_b1_placed['z'])**2
                    if dist_sq < 4:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # If all checks pass, add this B1 ball.
                b1_centers.append(p_b1)

    num_b1 = len(b1_centers)
    b1_value = num_b1 * 1

    total_value = b2_value + b1_value
    
    # Print the final equation with all its components.
    print(f"{num_b2} * 150 + {num_b1} * 1 = {b2_value} + {b1_value} = {total_value}")

solve_cutting_problem()
<<<2066>>>