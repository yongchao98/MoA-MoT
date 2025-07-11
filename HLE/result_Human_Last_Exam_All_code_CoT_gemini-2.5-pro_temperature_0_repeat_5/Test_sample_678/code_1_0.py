import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600},  # diameter 10m -> radius 5m
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000}
    }

    # --- Simulation Setup ---
    # Use a 1x1 meter grid for the room. We use a grid of points from (0,0) to (140, 110).
    total_points = (ROOM_W + 1) * (ROOM_H + 1)
    target_points = total_points * TARGET_COVERAGE_RATIO

    # Grid representing the room, True means covered
    coverage_grid = np.zeros((ROOM_H + 1, ROOM_W + 1), dtype=bool)

    # Possible center coordinates for scanners (multiples of 5)
    possible_locations = []
    for x in range(0, ROOM_W + 1, 5):
        for y in range(0, ROOM_H + 1, 5):
            possible_locations.append((x, y))

    # --- Pre-compute scanner masks for performance ---
    scanner_masks = {}
    y_coords, x_coords = np.ogrid[0:ROOM_H + 1, 0:ROOM_W + 1]

    print("Pre-computing scanner coverage masks...")
    for name, props in SCANNERS.items():
        for loc in possible_locations:
            px, py = loc
            if props['shape'] == 'circle':
                radius = props['radius']
                mask = (x_coords - px)**2 + (y_coords - py)**2 <= radius**2
            elif props['shape'] == 'square':
                side = props['side']
                half_side = side / 2
                mask = (np.abs(x_coords - px) <= half_side) & (np.abs(y_coords - py) <= half_side)
            scanner_masks[(name, loc)] = mask
    print("...done.")

    # --- Greedy Algorithm ---
    total_cost = 0
    placed_scanners = []
    equation_parts = []

    print("\nStarting optimization...")
    while np.sum(coverage_grid) < target_points:
        best_placement = None
        max_effectiveness = -1
        
        # Find the most cost-effective scanner to add
        for (name, loc), mask in scanner_masks.items():
            cost = SCANNERS[name]['cost']
            
            # Calculate how many *new* points this scanner would cover
            newly_covered_mask = mask & ~coverage_grid
            new_area = np.sum(newly_covered_mask)
            
            if new_area == 0:
                continue
            
            effectiveness = new_area / cost
            
            if effectiveness > max_effectiveness:
                max_effectiveness = effectiveness
                best_placement = (name, loc)

        if best_placement is None:
            print("No more effective placements found.")
            break

        # Add the best scanner found in this iteration
        name, loc = best_placement
        cost = SCANNERS[name]['cost']
        
        # Update coverage and cost
        mask_to_add = scanner_masks[best_placement]
        coverage_grid |= mask_to_add
        
        total_cost += cost
        placed_scanners.append(best_placement)
        equation_parts.append(str(cost))
        
        current_coverage_pts = np.sum(coverage_grid)
        coverage_percent = (current_coverage_pts / total_points) * 100
        
        print(f"Added {name} at {loc}. Total Cost: {total_cost}. Coverage: {coverage_percent:.2f}%")

    # --- Final Results ---
    print("\n--- Optimization Complete ---")

    num_c2 = sum(1 for s in placed_scanners if s[0] == 'C2')
    num_c1 = sum(1 for s in placed_scanners if s[0] == 'C1')
    num_r1 = sum(1 for s in placed_scanners if s[0] == 'R1')

    print(f"\nFinal Scanner Configuration:")
    print(f"- C2 (20m radius, 20000 cost): {num_c2} units")
    print(f"- C1 (5m radius, 1600 cost): {num_c1} units")
    print(f"- R1 (10m square, 2000 cost): {num_r1} units")

    final_coverage_pts = np.sum(coverage_grid)
    final_coverage_ratio = final_coverage_pts / total_points

    print(f"\nFinal Coverage: {final_coverage_ratio:.4f} (Target was {TARGET_COVERAGE_RATIO:.4f})")
    print(f"Total Area Points: {total_points}")
    print(f"Covered Area Points: {final_coverage_pts}")

    print("\nFinal Cost Calculation:")
    cost_equation = " + ".join(equation_parts)
    print(f"{cost_equation} = {total_cost}")

    print(f"\nOptimal Total Cost: {total_cost}")
    
    return total_cost

if __name__ == '__main__':
    solve_scanner_placement()