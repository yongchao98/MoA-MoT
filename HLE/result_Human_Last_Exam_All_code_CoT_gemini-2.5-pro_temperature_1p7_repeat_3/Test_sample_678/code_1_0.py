import numpy as np
import math

def solve_scanner_placement():
    """
    Calculates the optimal scanner placement using a greedy algorithm to minimize
    cost while achieving a target coverage level.
    """
    # --- Problem Definition ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    TARGET_COVERAGE_RATIO = 0.88

    TOTAL_POINTS = ROOM_W * ROOM_H
    TARGET_COVERAGE_POINTS = TOTAL_POINTS * TARGET_COVERAGE_RATIO

    # Scanner definitions
    scanners = {
        'C2': {'shape': 'circle', 'radius': 20, 'side': None, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'side': None, 'cost': 1600}, # diameter 10m -> radius 5m
        'R1': {'shape': 'square', 'side': 10, 'radius': None, 'cost': 2000}
    }

    # --- Grid Setup ---
    # Possible center locations for scanners (multiples of 5)
    possible_x = list(range(0, ROOM_W + 1, 5))
    possible_y = list(range(0, ROOM_H + 1, 5))
    placement_locations = [(x, y) for x in possible_x for y in possible_y]

    # Coverage grid (1 point per m^2) representing the room floor
    coverage_map = np.zeros((ROOM_H, ROOM_W), dtype=bool)

    # Coordinate arrays for efficient mask calculation
    gy, gx = np.ogrid[0:ROOM_H, 0:ROOM_W]

    # --- Helper Function for Mask Generation ---
    def get_scanner_mask(scanner_type_info, center_x, center_y):
        """Generates a boolean mask for a given scanner at a specific location."""
        if scanner_type_info['shape'] == 'circle':
            # Vectorized calculation for all points in a circle
            dist_sq = (gx - center_x)**2 + (gy - center_y)**2
            return dist_sq <= scanner_type_info['radius']**2
        elif scanner_type_info['shape'] == 'square':
            mask = np.zeros((ROOM_H, ROOM_W), dtype=bool)
            half_side = scanner_type_info['side'] / 2
            # Calculate boundaries, ensuring they are within the room
            x_min = max(0, round(center_x - half_side))
            x_max = min(ROOM_W, round(center_x + half_side))
            y_min = max(0, round(center_y - half_side))
            y_max = min(ROOM_H, round(center_y + half_side))
            mask[y_min:y_max, x_min:x_max] = True
            return mask

    # --- Greedy Algorithm ---
    total_cost = 0
    current_coverage_points = 0
    placed_scanners_list = []

    while current_coverage_points < TARGET_COVERAGE_POINTS:
        best_marginal_value = -1.0
        best_choice = None
        best_mask = None
        
        # Find the most cost-effective scanner to add in this step
        for s_type, s_info in scanners.items():
            cost = s_info['cost']
            for loc in placement_locations:
                center_x, center_y = loc
                
                scanner_mask = get_scanner_mask(s_info, center_x, center_y)
                
                # Find only the newly covered area by this scanner
                newly_covered_mask = np.logical_and(scanner_mask, np.logical_not(coverage_map))
                additional_coverage = np.sum(newly_covered_mask)
                
                if additional_coverage == 0:
                    continue
                
                marginal_value = additional_coverage / cost
                
                if marginal_value > best_marginal_value:
                    best_marginal_value = marginal_value
                    best_choice = {'type': s_type, 'cost': cost}
                    best_mask = scanner_mask

        if best_choice is None:
            print("Stopping: Could not find any scanner to add more coverage.")
            break
            
        # Add the chosen scanner to our solution
        placed_scanners_list.append(best_choice)
        total_cost += best_choice['cost']
        
        # Update the master coverage map by taking the union
        coverage_map = np.logical_or(coverage_map, best_mask)
        current_coverage_points = np.sum(coverage_map)

    # --- Final Output ---
    final_coverage_ratio = np.sum(coverage_map) / TOTAL_POINTS
    
    count_c2 = sum(1 for s in placed_scanners_list if s['type'] == 'C2')
    count_c1 = sum(1 for s in placed_scanners_list if s['type'] == 'C1')
    count_r1 = sum(1 for s in placed_scanners_list if s['type'] == 'R1')

    print("--- Optimization Finished ---")
    print(f"Final Coverage Achieved: {final_coverage_ratio:.2%}")
    print(f"\nOptimal Scanner Configuration:")
    print(f"- C2 Scanners (20m radius, 20000 cost): {count_c2}")
    print(f"- C1 Scanners (5m radius, 1600 cost):  {count_c1}")
    print(f"- R1 Scanners (10m square, 2000 cost): {count_r1}")
    
    # Final equation as requested
    print(f"\nFinal Equation for Total Cost:")
    print(f"Total Cost = ({count_c2} * {scanners['C2']['cost']}) + ({count_c1} * {scanners['C1']['cost']}) + ({count_r1} * {scanners['R1']['cost']}) = {total_cost}")

    print(f"\nOptimal Total Cost:")
    print(total_cost)

    return total_cost

if __name__ == '__main__':
    solve_scanner_placement()
