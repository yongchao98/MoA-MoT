import numpy as np
import math

def solve_museum_coverage():
    """
    Solves the museum scanner placement problem to achieve at least 88% coverage
    with the minimum cost using a greedy algorithm.
    """
    # 1. Define constants and room parameters
    ROOM_W = 140  # Room width in meters
    ROOM_H = 110  # Room height in meters
    GRID_MULTIPLE = 5
    TARGET_COVERAGE_RATIO = 0.88

    TOTAL_AREA = ROOM_W * ROOM_H
    TARGET_AREA = TOTAL_AREA * TARGET_COVERAGE_RATIO

    print(f"Room Dimensions: {ROOM_W}m x {ROOM_H}m (Total Area: {TOTAL_AREA} m^2)")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.0%} (Target Area: {TARGET_AREA} m^2)\n")


    # 2. Define scanner properties
    scanners = {
        'C2': {'type': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'type': 'circle', 'radius': 5, 'cost': 1600}, # 10m diameter
        'R1': {'type': 'square', 'side': 10, 'cost': 2000}
    }

    # 3. Model the room and potential placements
    # Create a grid where each cell is 1x1m, representing the floor
    # We use (height, width) for numpy array shape, which corresponds to (y, x)
    coverage_map = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    
    # Generate a list of all possible center coordinates for scanners
    possible_centers = [
        (x, y) for x in range(0, ROOM_W + 1, GRID_MULTIPLE)
        for y in range(0, ROOM_H + 1, GRID_MULTIPLE)
    ]

    # Pre-compute all possible scanner coverage masks to speed up the process
    # This avoids recalculating the shape of each scanner at every step
    print("Pre-computing all possible scanner placement masks...")
    all_masks = {}
    grid_y, grid_x = np.ogrid[:ROOM_H, :ROOM_W] # Coordinate grid
    for name, info in scanners.items():
        for cx, cy in possible_centers:
            if info['type'] == 'circle':
                r = info['radius']
                mask = (grid_x - cx)**2 + (grid_y - cy)**2 <= r**2
            elif info['type'] == 'square':
                s_half = info['side'] / 2
                mask = (np.abs(grid_x - cx) <= s_half) & (np.abs(grid_y - cy) <= s_half)
            all_masks[(name, cx, cy)] = mask
    print("Pre-computation complete.\n")


    # 4. Greedy algorithm to place scanners
    total_cost = 0
    placed_scanners = []
    current_coverage = 0
    step = 1

    print("Starting optimization process...\n")
    while current_coverage < TARGET_AREA:
        best_value = -1
        best_choice = None
        best_new_area = 0

        # Find the best scanner to place in this step
        for (name, cx, cy), mask in all_masks.items():
            cost = scanners[name]['cost']
            
            # Calculate the new area this scanner would cover
            # by looking at the part of its mask that doesn't overlap with existing coverage
            newly_covered_mask = mask & (~coverage_map)
            new_area = np.sum(newly_covered_mask)

            if new_area > 0:
                # Value is new area per unit cost (cost-effectiveness)
                value = new_area / cost
                if value > best_value:
                    best_value = value
                    best_choice = (name, cx, cy)
                    best_new_area = new_area
        
        if best_choice is None:
            print("No further coverage improvement is possible.")
            break

        # Place the best scanner found
        name, cx, cy = best_choice
        cost = scanners[name]['cost']
        mask_to_add = all_masks[best_choice]

        # Update state
        coverage_map |= mask_to_add
        total_cost += cost
        current_coverage = np.sum(coverage_map)
        placed_scanners.append({'name': name, 'center': (cx, cy), 'cost': cost})
        
        print(f"Step {step}: Placed {name} scanner at ({cx}, {cy}). "
              f"Added {best_new_area} m^2. "
              f"Cost: {total_cost}. "
              f"Coverage: {current_coverage / TOTAL_AREA:.2%}")
        step += 1

    # 5. Report the final results
    print("\n--- Optimization Complete ---\n")

    n_c2 = sum(1 for s in placed_scanners if s['name'] == 'C2')
    n_c1 = sum(1 for s in placed_scanners if s['name'] == 'C1')
    n_r1 = sum(1 for s in placed_scanners if s['name'] == 'R1')

    print("Final Scanner Configuration:")
    print(f"- {n_c2} C2 scanners (Circle, 20m radius)")
    print(f"- {n_c1} C1 scanners (Circle, 5m radius)")
    print(f"- {n_r1} R1 scanners (Square, 10m side)")
    
    final_coverage = np.sum(coverage_map)
    final_coverage_ratio = final_coverage / TOTAL_AREA
    print(f"\nFinal Coverage: {final_coverage} / {TOTAL_AREA} m^2 = {final_coverage_ratio:.2%}")
    
    print("\nFinal Cost Calculation:")
    cost_c2 = scanners['C2']['cost']
    cost_c1 = scanners['C1']['cost']
    cost_r1 = scanners['R1']['cost']
    
    # Print the equation with each number explicitly
    print(f"Total Cost = ({n_c2} * {cost_c2}) + ({n_c1} * {cost_c1}) + ({n_r1} * {cost_r1})")
    print(f"Total Cost = {n_c2 * cost_c2} + {n_c1 * cost_c1} + {n_r1 * cost_r1}")
    print(f"Total Cost = {total_cost}")

    print(f"\n<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_museum_coverage()