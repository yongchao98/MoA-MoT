import numpy as np
import math
from collections import Counter

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Setup
    ROOM_W = 140
    ROOM_H = 110
    TOTAL_CELLS = ROOM_W * ROOM_H
    TARGET_COVERAGE_RATIO = 0.88
    TARGET_COVERED_CELLS = TOTAL_CELLS * TARGET_COVERAGE_RATIO

    scanners_spec = {
        'C2': {'type': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'type': 'circle', 'radius': 5, 'cost': 1600},
        'R1': {'type': 'square', 'side': 10, 'cost': 2000},
    }

    # Create a grid of cell-center coordinates for coverage calculation
    y_coords, x_coords = np.mgrid[0.5:ROOM_H, 0.5:ROOM_W]

    # 2. Heuristic Base Configuration (5 C2 scanners)
    # This is a strong starting point based on preliminary analysis
    placements = [
        ('C2', 25, 25), ('C2', 25, 85),
        ('C2', 115, 25), ('C2', 115, 85),
        ('C2', 70, 55)
    ]

    # Calculate initial coverage mask and cost
    covered_mask = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    current_cost = 0
    for name, cx, cy in placements:
        spec = scanners_spec[name]
        current_cost += spec['cost']
        if spec['type'] == 'circle':
            dist_sq = (x_coords - cx)**2 + (y_coords - cy)**2
            covered_mask |= (dist_sq <= spec['radius']**2)
    
    covered_cells = np.sum(covered_mask)

    # 3. Greedy Loop to add fillers
    possible_locations = []
    for x in range(0, ROOM_W + 1, 5):
        for y in range(0, ROOM_H + 1, 5):
            possible_locations.append((x, y))

    addon_scanner_names = ['R1', 'C1']

    while covered_cells < TARGET_COVERED_CELLS:
        best_placement = None
        max_bang_for_buck = -1
        best_new_mask = None

        for name in addon_scanner_names:
            spec = scanners_spec[name]
            cost = spec['cost']
            
            for cx, cy in possible_locations:
                # Avoid placing on top of existing scanners
                if any(p[1] == cx and p[2] == cy for p in placements):
                    continue

                # Calculate mask for this potential scanner
                if spec['type'] == 'circle':
                    dist_sq = (x_coords - cx)**2 + (y_coords - cy)**2
                    scanner_mask = (dist_sq <= spec['radius']**2)
                else: # square
                    half_side = spec['side'] / 2.0
                    scanner_mask = (np.abs(x_coords - cx) <= half_side) & \
                                   (np.abs(y_coords - cy) <= half_side)

                # Calculate newly covered area and bang-for-buck
                newly_covered_mask = scanner_mask & ~covered_mask
                new_cells = np.sum(newly_covered_mask)
                
                if new_cells == 0:
                    continue

                bang_for_buck = new_cells / cost

                if bang_for_buck > max_bang_for_buck:
                    max_bang_for_buck = bang_for_buck
                    best_placement = (name, cx, cy)
                    best_new_mask = newly_covered_mask

        if best_placement is None:
            print("Warning: Could not find any more valuable placements. Stopping.")
            break
            
        # Add the best found scanner
        placements.append(best_placement)
        spec = scanners_spec[best_placement[0]]
        current_cost += spec['cost']
        covered_mask |= best_new_mask
        covered_cells = np.sum(covered_mask)

    # 4. Output Results
    print("Optimal scanner configuration found:")
    
    counts = Counter(p[0] for p in placements)
    
    print(f"\n- C2 Scanners (Cost {scanners_spec['C2']['cost']}): {counts['C2']}")
    print(f"- C1 Scanners (Cost {scanners_spec['C1']['cost']}): {counts['C1']}")
    print(f"- R1 Scanners (Cost {scanners_spec['R1']['cost']}): {counts['R1']}")
    
    final_coverage = covered_cells / TOTAL_CELLS
    print(f"\nFinal Coverage: {final_coverage:.2%} (Target: {TARGET_COVERAGE_RATIO:.2%})")
    
    cost_parts = []
    total_cost_calc = 0
    
    if counts['C2'] > 0:
        cost = counts['C2'] * scanners_spec['C2']['cost']
        cost_parts.append(f"{counts['C2']} * {scanners_spec['C2']['cost']}")
        total_cost_calc += cost
    if counts['C1'] > 0:
        cost = counts['C1'] * scanners_spec['C1']['cost']
        cost_parts.append(f"{counts['C1']} * {scanners_spec['C1']['cost']}")
        total_cost_calc += cost
    if counts['R1'] > 0:
        cost = counts['R1'] * scanners_spec['R1']['cost']
        cost_parts.append(f"{counts['R1']} * {scanners_spec['R1']['cost']}")
        total_cost_calc += cost
        
    equation = " + ".join(cost_parts)
    print(f"\nTotal Cost Calculation: {equation} = {total_cost_calc}")
    print("\nOptimal Total Cost:")
    print(f"<<<{total_cost_calc}>>>")

solve_scanner_placement()