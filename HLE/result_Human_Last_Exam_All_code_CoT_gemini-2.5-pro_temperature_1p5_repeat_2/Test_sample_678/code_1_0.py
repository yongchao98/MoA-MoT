import numpy as np
import math
from collections import Counter

def solve_scanner_placement():
    """
    Solves the scanner placement optimization problem using a greedy heuristic.
    """
    # --- Configuration ---
    ROOM_W = 140
    ROOM_H = 110
    COVERAGE_TARGET_RATIO = 0.88
    GRID_STEP = 5

    ROOM_AREA = ROOM_W * ROOM_H
    TARGET_COVERAGE_AREA = ROOM_AREA * COVERAGE_TARGET_RATIO

    # --- Scanner Definitions ---
    scanners = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10,  'cost': 2000}
    }

    # 1. Generate potential center locations
    possible_centers = []
    for x in range(0, ROOM_W + 1, GRID_STEP):
        for y in range(0, ROOM_H + 1, GRID_STEP):
            possible_centers.append((x, y))

    # 2. Create the coverage grid for the room (1m resolution)
    coverage_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    
    # Pre-generate coordinates for mask creation, this is more efficient
    y_coords, x_coords = np.ogrid[:ROOM_H, :ROOM_W]

    # --- Greedy Algorithm ---
    total_cost = 0
    placed_scanners = []
    current_coverage_area = 0

    print("Running optimization, this may take a moment...")

    while current_coverage_area < TARGET_COVERAGE_AREA:
        best_value = -1
        best_placement_info = None
        
        # Find the best scanner to add in this iteration
        for center in possible_centers:
            cx, cy = center
            for s_type, s_info in scanners.items():
                
                # Generate the mask for this potential placement on the fly
                if s_info['shape'] == 'circle':
                    radius = s_info['radius']
                    dist_sq = (x_coords - cx)**2 + (y_coords - cy)**2
                    potential_mask = dist_sq <= radius**2
                elif s_info['shape'] == 'square':
                    side = s_info['side']
                    # Covers a 10x10 area of cells
                    mask_x = (x_coords >= cx - side / 2) & (x_coords < cx + side / 2)
                    mask_y = (y_coords >= cy - side / 2) & (y_coords < cy + side / 2)
                    potential_mask = mask_x & mask_y
                
                # Calculate new coverage this scanner would add
                new_coverage_mask = np.logical_and(potential_mask, ~coverage_grid)
                newly_covered_area = np.sum(new_coverage_mask)
                
                if newly_covered_area == 0:
                    continue
                
                # Calculate the value: new area per unit cost
                value = newly_covered_area / s_info['cost']
                
                if value > best_value:
                    best_value = value
                    best_placement_info = {
                        'type': s_type,
                        'center': center,
                        'cost': s_info['cost'],
                        'mask': potential_mask
                    }
        
        if best_placement_info is None:
            print("Warning: No more coverage can be added, but target not met.")
            break

        # Add the best scanner found in the iteration
        placed_scanners.append(best_placement_info)
        total_cost += best_placement_info['cost']
        
        # Update the main coverage grid
        coverage_grid = np.logical_or(coverage_grid, best_placement_info['mask'])
        
        # Update the total covered area
        current_coverage_area = np.sum(coverage_grid)

    # --- Final Results ---
    print("\n--- Optimization Complete ---")
    print(f"Target coverage: {COVERAGE_TARGET_RATIO:.2%} ({int(TARGET_COVERAGE_AREA)} m^2)")
    print(f"Achieved coverage: {current_coverage_area / ROOM_AREA:.2%} ({int(current_coverage_area)} m^2)")
    
    scanner_counts = Counter(s['type'] for s in placed_scanners)
    
    print("\nOptimal Scanner Configuration:")
    cost_equation_parts = []
    for s_type, count in sorted(scanner_counts.items()):
        cost = scanners[s_type]['cost']
        sub_total = count * cost
        cost_str = f"{count} * {cost}"
        print(f"- {count} x {s_type} scanners: {cost_str} = {sub_total}")
        cost_equation_parts.append(str(sub_total))

    total_cost_str = " + ".join(cost_equation_parts)
    print(f"\nTotal Cost = {total_cost_str} = {total_cost}")
    
    print(f"\nFinal Answer:")
    print(f"<<<{total_cost}>>>")


solve_scanner_placement()
