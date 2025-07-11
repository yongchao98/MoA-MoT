import numpy as np

def solve_scanner_placement():
    """
    Solves the museum scanner placement optimization problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    PLACEMENT_STEP = 5 # meters
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = [
        {'name': 'C2', 'cost': 20000, 'type': 'circle', 'param': 20}, # param is radius
        {'name': 'R1', 'cost': 2000, 'type': 'square', 'param': 10},  # param is side
        {'name': 'C1', 'cost': 1600, 'type': 'circle', 'param': 5},   # param is radius (from 10m diameter)
    ]
    
    # Sort scanners by a simple cost-effectiveness metric to potentially speed up finding the best move
    # This is a minor heuristic and doesn't change the core greedy logic
    for s in SCANNERS:
        if s['type'] == 'circle':
            area = np.pi * s['param']**2
        else: # square
            area = s['param']**2
        s['base_effectiveness'] = s['cost'] / area
    SCANNERS.sort(key=lambda x: x['base_effectiveness'])


    # 1. Setup grids and constants
    total_pixels = ROOM_W * ROOM_H
    target_pixels_to_cover = total_pixels * TARGET_COVERAGE_RATIO
    
    # Represents the room as a grid of 1x1m squares
    coverage_map = np.zeros((ROOM_H, ROOM_W), dtype=bool)

    # Generate possible placement locations
    placement_coords = []
    for x in range(0, ROOM_W + 1, PLACEMENT_STEP):
        for y in range(0, ROOM_H + 1, PLACEMENT_STEP):
            placement_coords.append((x, y))

    # Pre-calculate coverage masks for all possible placements to speed up the main loop
    precalculated_masks = {}
    yy, xx = np.mgrid[0.5:ROOM_H, 0.5:ROOM_W]
    for sc in SCANNERS:
        for (cx, cy) in placement_coords:
            if sc['type'] == 'circle':
                radius = sc['param']
                # Equation of a circle: (x-cx)^2 + (y-cy)^2 <= r^2
                mask = (xx - cx)**2 + (yy - cy)**2 <= radius**2
            elif sc['type'] == 'square':
                side = sc['param']
                x_min, x_max = round(cx - side / 2), round(cx + side / 2)
                y_min, y_max = round(cy - side / 2), round(cy + side / 2)
                # Create a boolean mask for the square
                mask = (xx >= x_min) & (xx < x_max) & (yy >= y_min) & (yy < y_max)
            
            # Clip the mask to the room's boundaries
            bounded_mask = np.zeros_like(coverage_map, dtype=bool)
            bounded_mask = mask
            precalculated_masks[(sc['name'], cx, cy)] = bounded_mask

    # 2. Greedy search loop
    total_cost = 0
    total_covered_pixels = 0
    solution = []
    used_locations = set()

    print("--- Starting Optimization ---")
    while total_covered_pixels < target_pixels_to_cover:
        best_move = None
        min_cost_per_pixel = float('inf')

        for (cx, cy) in placement_coords:
            if (cx, cy) in used_locations:
                continue

            for sc in SCANNERS:
                potential_mask = precalculated_masks[(sc['name'], cx, cy)]
                
                # Calculate new coverage by finding pixels in the potential mask 
                # that are not already covered in the main coverage_map.
                newly_covered_mask = potential_mask & ~coverage_map
                new_pixels = np.sum(newly_covered_mask)
                
                if new_pixels == 0:
                    continue

                cost_per_pixel = sc['cost'] / new_pixels
                
                if cost_per_pixel < min_cost_per_pixel:
                    min_cost_per_pixel = cost_per_pixel
                    best_move = {
                        'scanner': sc,
                        'pos': (cx, cy),
                        'new_mask': newly_covered_mask,
                        'new_pixels': new_pixels
                    }

        if best_move is None:
            print("No further coverage improvement possible.")
            break
            
        # Apply the best move found in the iteration
        scanner_info = best_move['scanner']
        pos = best_move['pos']
        
        # Update state
        total_cost += scanner_info['cost']
        coverage_map |= best_move['new_mask'] 
        total_covered_pixels = np.sum(coverage_map)
        used_locations.add(pos)
        solution.append({'name': scanner_info['name'], 'pos': pos, 'cost': scanner_info['cost']})
        
        # Print progress
        print(f"Placed {scanner_info['name']} at {pos}. Cost: {scanner_info['cost']:5d}. "
              f"Coverage: {total_covered_pixels / total_pixels:.2%}. "
              f"Total Cost: {total_cost}")

    # 4. Final results
    print("\n--- Optimization Complete ---")
    
    counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for item in solution:
        counts[item['name']] += 1

    print("Optimal Scanner Configuration:")
    final_cost_strings = []
    for name in ['C2', 'C1', 'R1']:
        count = counts[name]
        if count > 0:
            cost = next(sc['cost'] for sc in SCANNERS if sc['name'] == name)
            print(f"- Type {name}: {count} units")
            final_cost_strings.append(f"{count} * {cost}")
    
    final_coverage_ratio = total_covered_pixels / total_pixels
    print(f"\nFinal Coverage: {total_covered_pixels} / {total_pixels} pixels = {final_coverage_ratio:.4f}")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.4f} (Reached: {final_coverage_ratio >= TARGET_COVERAGE_RATIO})")

    print("\nFinal Cost Calculation:")
    print(f"{' + '.join(final_cost_strings)} = {total_cost}")

    print(f"\n<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()