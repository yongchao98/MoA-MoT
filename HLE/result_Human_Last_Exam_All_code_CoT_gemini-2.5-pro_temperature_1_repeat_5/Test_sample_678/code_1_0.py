import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Step 1: Configuration & Initialization ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    PLACEMENT_STEP = 5
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'side': None, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'side': None, 'cost': 1600},
        'R1': {'shape': 'square', 'side': 10, 'radius': None, 'cost': 2000}
    }

    TOTAL_AREA_POINTS = ROOM_WIDTH * ROOM_HEIGHT
    TARGET_COVERAGE_POINTS = math.ceil(TOTAL_AREA_POINTS * TARGET_COVERAGE_RATIO)

    # A grid representing the room, where each cell is a 1x1m square.
    coverage_grid = np.zeros((ROOM_HEIGHT, ROOM_WIDTH), dtype=bool)

    # --- Step 2: Identify Possible Placements ---
    possible_centers = set()
    for x in range(0, ROOM_WIDTH + 1, PLACEMENT_STEP):
        for y in range(0, ROOM_HEIGHT + 1, PLACEMENT_STEP):
            possible_centers.add((x, y))

    # --- Step 3: Pre-calculate Scanner Footprints for speed ---
    footprints = {}
    room_y_coords, room_x_coords = np.mgrid[0:ROOM_HEIGHT, 0:ROOM_WIDTH]

    for center_x, center_y in possible_centers:
        for name, props in SCANNERS.items():
            if props['shape'] == 'circle':
                radius = props['radius']
                # Equation of a circle: (x-h)^2 + (y-k)^2 <= r^2
                dist_sq = (room_x_coords - center_x)**2 + (room_y_coords - center_y)**2
                mask = dist_sq <= radius**2
            elif props['shape'] == 'square':
                side = props['side']
                half_side = side / 2.0
                mask = (room_x_coords >= center_x - half_side) & \
                       (room_x_coords < center_x + half_side) & \
                       (room_y_coords >= center_y - half_side) & \
                       (room_y_coords < center_y + half_side)
            footprints[(name, (center_x, center_y))] = mask
    
    # --- Step 4: Greedy Algorithm Loop ---
    total_cost = 0
    placed_scanners_count = {'C2': 0, 'C1': 0, 'R1': 0}
    
    while np.sum(coverage_grid) < TARGET_COVERAGE_POINTS:
        current_coverage_pts = np.sum(coverage_grid)
        
        best_option = {'value': -1}

        # Iterate over centers that are not yet used
        for center in list(possible_centers):
            # For each center, evaluate the 3 scanner types
            for name, props in SCANNERS.items():
                footprint_mask = footprints[(name, center)]
                
                # Calculate new coverage by finding points in footprint NOT in current coverage_grid
                newly_covered_points = np.sum(footprint_mask & ~coverage_grid)

                if newly_covered_points == 0:
                    continue
                
                cost = props['cost']
                value = newly_covered_points / cost

                if value > best_option['value']:
                    best_option = {
                        'value': value,
                        'name': name,
                        'center': center,
                        'cost': cost,
                    }
        
        if best_option['value'] <= 0:
            print("\nStopping: No cost-effective placement adds new coverage.")
            break

        # Place the best scanner found in this iteration
        best_name = best_option['name']
        best_center = best_option['center']
        best_cost = best_option['cost']
        
        # Add to totals
        total_cost += best_cost
        placed_scanners_count[best_name] += 1
        
        # Update coverage grid
        best_footprint = footprints[(best_name, best_center)]
        coverage_grid |= best_footprint
        
        # This center location is now used and cannot be used again
        possible_centers.remove(best_center)

    # --- Step 5: Final Result ---
    print("Optimization Complete.")
    final_coverage = np.sum(coverage_grid)
    final_coverage_ratio = final_coverage / TOTAL_AREA_POINTS
    
    print(f"Final Coverage Achieved: {final_coverage_ratio:.4f} (Target was {TARGET_COVERAGE_RATIO:.2f})")
    
    print("\nFinal Cost Calculation:")
    num_c2 = placed_scanners_count['C2']
    num_c1 = placed_scanners_count['C1']
    num_r1 = placed_scanners_count['R1']
    
    # Print the final equation as requested
    print(f"{num_c2} * 20000 + {num_c1} * 1600 + {num_r1} * 2000 = {total_cost}")

if __name__ == '__main__':
    solve_scanner_placement()