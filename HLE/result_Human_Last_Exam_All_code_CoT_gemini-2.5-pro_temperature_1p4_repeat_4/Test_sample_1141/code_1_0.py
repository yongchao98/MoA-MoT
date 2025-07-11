import math

def solve_cutting_problem():
    """
    Solves the billet cutting optimization problem based on the user's formulation.
    """
    
    # --- Problem Definition ---
    # Prices and dimensions (in 0.5 cm grid units)
    pieces_info = {
        'B2': {'price': 150, 'radius': 4},
        'B1': {'price': 1, 'radius': 1},
        'T1': {'price': 5} 
    }
    
    # Center coordinate ranges for T1 and B1
    b1_bounds = {'x': (1, 31), 'y': (1, 21), 'z': (1, 7)}

    # --- Constraint Checking Functions (as per user formulation) ---
    def dist_sq(p1, p2):
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def check_overlap(item1_type, item1_pos, item2_type, item2_pos):
        # Ball to Ball
        if item1_type in ['B1', 'B2'] and item2_type in ['B1', 'B2']:
            r1 = pieces_info[item1_type]['radius']
            r2 = pieces_info[item2_type]['radius']
            min_dist_sq = (r1 + r2)**2
            return dist_sq(item1_pos, item2_pos) < min_dist_sq
        # All other constraints involve T1, which we concluded are not in the optimal solution.
        # This function is simplified to only handle the relevant cases for the chosen strategy.
        return False

    # --- Strategy: Maximize B2s, then fill with B1s ---
    
    # 1. Place B2 balls in a grid
    # This grid arrangement maximizes the number of B2s in the given billet.
    placed_b2 = []
    for x in [4, 12, 20, 28]:
        for y in [4, 12]:
            placed_b2.append(('B2', (x, y, 4)))

    # 2. Greedily add B1 balls into the remaining space
    placed_b1 = []
    
    # Iterate through all possible B1 center coordinates in a specific order
    for z in range(b1_bounds['z'][0], b1_bounds['z'][1] + 1):
        for y in range(b1_bounds['y'][0], b1_bounds['y'][1] + 1):
            for x in range(b1_bounds['x'][0], b1_bounds['x'][1] + 1):
                new_b1_pos = (x, y, z)
                is_valid = True
                
                # Check for overlap with the placed B2 balls
                for b2_type, b2_pos in placed_b2:
                    if check_overlap('B1', new_b1_pos, b2_type, b2_pos):
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # Check for overlap with already placed B1 balls
                for b1_type, b1_pos in placed_b1:
                    if check_overlap('B1', new_b1_pos, b1_type, b1_pos):
                        is_valid = False
                        break
                
                # If no overlaps, place this B1 ball
                if is_valid:
                    placed_b1.append(('B1', new_b1_pos))

    # --- Calculate Final Result ---
    num_b2 = len(placed_b2)
    price_b2 = pieces_info['B2']['price']
    
    num_b1 = len(placed_b1)
    price_b1 = pieces_info['B1']['price']
    
    num_t1 = 0 # As reasoned, T1s are not included
    price_t1 = pieces_info['T1']['price']
    
    total_value = (num_b2 * price_b2) + (num_b1 * price_b1) + (num_t1 * price_t1)
    
    # --- Print Final Equation ---
    print(f"The highest value solution found is by cutting {num_b2} B2 pieces, {num_b1} B1 pieces, and {num_t1} T1 pieces.")
    print("The calculation for the total value is:")
    print(f"{num_b2} * {price_b2} + {num_b1} * {price_b1} + {num_t1} * {price_t1} = {total_value}")


solve_cutting_problem()
<<<1680>>>