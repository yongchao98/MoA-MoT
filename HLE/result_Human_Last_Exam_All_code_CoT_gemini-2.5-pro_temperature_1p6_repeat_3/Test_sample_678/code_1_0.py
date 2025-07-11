import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_W, ROOM_H = 140, 110
    ROOM_AREA = ROOM_W * ROOM_H
    TARGET_COVERAGE_RATIO = 0.88
    TARGET_AREA = ROOM_AREA * TARGET_COVERAGE_RATIO
    GRID_STEP = 5

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600},
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000}
    }

    # --- Pre-computation: Create boolean masks for each scanner's coverage area ---
    scanner_masks = {}
    for name, props in SCANNERS.items():
        if props['shape'] == 'circle':
            r = props['radius']
            size = 2 * r + 1
            mask = np.zeros((size, size), dtype=bool)
            center = r
            for y in range(size):
                for x in range(size):
                    if (x - center)**2 + (y - center)**2 <= r**2:
                        mask[y, x] = True
            scanner_masks[name] = mask
        elif props['shape'] == 'square':
            s_half = props['side'] // 2
            size = 2 * s_half + 1
            mask = np.ones((size, size), dtype=bool)
            scanner_masks[name] = mask

    # --- Main Logic ---
    
    # 1. Initialization
    # A grid representing each 1x1m square in the room
    room_grid = np.zeros((ROOM_H + 1, ROOM_W + 1), dtype=bool)
    
    # Generate all possible scanner center locations
    possible_x = range(0, ROOM_W + 1, GRID_STEP)
    possible_y = range(0, ROOM_H + 1, GRID_STEP)
    possible_centers = [(x, y) for x in possible_x for y in possible_y]
    
    total_cost = 0
    covered_area = 0
    placements = []
    
    print(f"Room Area: {ROOM_AREA} m^2")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO * 100:.2f}% ({TARGET_AREA:.0f} m^2)")
    print("-" * 40)

    step = 1
    # 2. Greedy Loop: Continue adding scanners until target coverage is met
    while covered_area < TARGET_AREA:
        best_scanner_choice = None
        max_effectiveness = -1.0

        # Evaluate every scanner type at every possible location
        for name, props in SCANNERS.items():
            cost = props['cost']
            mask = scanner_masks[name]
            h, w = mask.shape
            offset_y, offset_x = h // 2, w // 2

            for cx, cy in possible_centers:
                # This loop calculates how many *new* cells a scanner would cover
                new_cells_count = 0
                for dy_mask in range(h):
                    for dx_mask in range(w):
                        # Check if the point is within the scanner's shape
                        if mask[dy_mask, dx_mask]:
                            # Find the corresponding absolute coordinates in the room
                            px, py = (cx - offset_x + dx_mask), (cy - offset_y + dy_mask)
                            # Check if the point is inside the room and not already covered
                            if (0 <= px <= ROOM_W and 0 <= py <= ROOM_H and not room_grid[py, px]):
                                new_cells_count += 1
                
                if new_cells_count == 0:
                    continue
                
                effectiveness = new_cells_count / cost

                if effectiveness > max_effectiveness:
                    max_effectiveness = effectiveness
                    best_scanner_choice = {
                        'name': name,
                        'center': (cx, cy),
                        'cost': cost,
                        'new_cells': new_cells_count
                    }
        
        if best_scanner_choice is None:
            print("No more effective placements found, stopping.")
            break

        # 3. Add the best scanner found in this iteration to the solution
        name = best_scanner_choice['name']
        center = best_scanner_choice['center']
        cost = best_scanner_choice['cost']
        new_cells = best_scanner_choice['new_cells']

        # Update totals
        total_cost += cost
        covered_area += new_cells
        placements.append(best_scanner_choice)
        
        # 4. Update the room grid to mark the newly covered area
        cx, cy = center
        mask = scanner_masks[name]
        h, w = mask.shape
        offset_y, offset_x = h // 2, w // 2
        for dy_mask in range(h):
            for dx_mask in range(w):
                if mask[dy_mask, dx_mask]:
                    px, py = (cx - offset_x + dx_mask), (cy - offset_y + dy_mask)
                    if 0 <= px <= ROOM_W and 0 <= py <= ROOM_H:
                        room_grid[py, px] = True
        
        coverage_percent = (covered_area / ROOM_AREA) * 100
        print(f"Step {step}: Adding {name} at {center} (Cost: {cost})")
        print(f"  > Added {new_cells} m^2 of new coverage.")
        print(f"  > Total Cost: {total_cost}, Total Coverage: {covered_area} / {TARGET_AREA:.0f} m^2 ({coverage_percent:.2f}%)")
        print("-" * 40)
        step += 1
        
    # --- Final Output ---
    print("\n--- Optimization Complete ---")
    final_coverage_percent = (covered_area / ROOM_AREA) * 100
    print(f"Final solution achieved coverage of {covered_area} m^2 ({final_coverage_percent:.2f}%), meeting the target of {TARGET_AREA:.0f} m^2.")
    
    counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for p in placements:
        counts[p['name']] += 1
        
    print("\nFinal Scanner Configuration:")
    print(f"  - C2 Scanners (Cost {SCANNERS['C2']['cost']} each): {counts['C2']}")
    print(f"  - C1 Scanners (Cost {SCANNERS['C1']['cost']} each): {counts['C1']}")
    print(f"  - R1 Scanners (Cost {SCANNERS['R1']['cost']} each): {counts['R1']}")
    
    print("\nFinal cost equation:")
    print(f"({counts['C2']} * {SCANNERS['C2']['cost']}) + ({counts['C1']} * {SCANNERS['C1']['cost']}) + ({counts['R1']} * {SCANNERS['R1']['cost']}) = {total_cost}")
    
    print(f"\nThe optimal total cost is: {total_cost}")

# Execute the solver
if __name__ == '__main__':
    solve_scanner_placement()
