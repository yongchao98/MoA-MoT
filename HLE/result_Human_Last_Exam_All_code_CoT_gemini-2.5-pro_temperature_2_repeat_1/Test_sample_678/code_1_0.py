import numpy as np
import math
import time

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Constants ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    TARGET_COVERAGE_RATIO = 0.88
    
    # Scanner properties: Type -> {shape, dimension, cost}
    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10,   'cost': 2000}
    }
    
    # --- Initialization ---
    total_area = ROOM_WIDTH * ROOM_HEIGHT
    target_covered_area = total_area * TARGET_COVERAGE_RATIO
    
    print(f"Room Area: {total_area} sq.m.")
    print(f"Target Coverage (>= {TARGET_COVERAGE_RATIO * 100}%): {target_covered_area:.2f} sq.m.\n")
    
    # The coverage_map represents the room. 0=uncovered, 1=covered.
    coverage_map = np.zeros((ROOM_HEIGHT, ROOM_WIDTH), dtype=bool)
    
    placed_scanners = []
    total_cost = 0
    iteration = 0

    # Create a list of all possible center coordinates for the scanners
    placement_coords = []
    for x in range(0, ROOM_WIDTH + 1, 5):
        for y in range(0, ROOM_HEIGHT + 1, 5):
            placement_coords.append((x, y))

    # Pre-calculate meshgrid for distance calculations
    yy, xx = np.mgrid[:ROOM_HEIGHT, :ROOM_WIDTH]
    
    start_time = time.time()
    # --- Greedy Algorithm Loop ---
    while np.sum(coverage_map) < target_covered_area:
        iteration += 1
        best_choice = {'value': -1}

        # Find the best scanner to add in this step
        for scanner_type, props in SCANNERS.items():
            for cx, cy in placement_coords:
                # Generate a boolean mask for the current scanner's coverage area
                if props['shape'] == 'circle':
                    mask = (xx - cx + 0.5)**2 + (yy - cy + 0.5)**2 <= props['radius']**2
                else: # square
                    half_side = props['side'] / 2
                    mask = (np.abs(xx - cx + 0.5) < half_side) & (np.abs(yy - cy + 0.5) < half_side)
                
                # Calculate the newly covered area by this scanner
                # This is the area the scanner covers (`mask`) that is not already covered (`~coverage_map`)
                additional_coverage = np.sum(mask & ~coverage_map)

                if additional_coverage > 0:
                    cost = props['cost']
                    value = additional_coverage / cost
                    if value > best_choice['value']:
                        best_choice.update({
                            'value': value,
                            'type': scanner_type,
                            'location': (cx, cy),
                            'mask': mask,
                            'cost': cost
                        })

        if best_choice['value'] == -1:
            print("No further coverage is possible. Stopping.")
            break

        # Add the best scanner from this iteration to the plan
        scanner_type = best_choice['type']
        location = best_choice['location']
        cost = best_choice['cost']

        total_cost += cost
        coverage_map |= best_choice['mask'] # Update the map using bitwise OR
        current_coverage_area = np.sum(coverage_map)
        placed_scanners.append({'type': scanner_type, 'cost': cost})

        print(f"Step {iteration}: Added {scanner_type} at {location}. "
              f"Coverage: {current_coverage_area/total_area:.4f}, "
              f"Cost: {total_cost}")

    end_time = time.time()
    
    # --- Final Results ---
    print("\n" + "="*40)
    print("Optimization Complete")
    print(f"Total execution time: {end_time - start_time:.2f} seconds")
    print("="*40)

    scanner_counts = {st: 0 for st in SCANNERS.keys()}
    for scanner in placed_scanners:
        scanner_counts[scanner['type']] += 1

    n_c2 = scanner_counts['C2']
    n_c1 = scanner_counts['C1']
    n_r1 = scanner_counts['R1']
    
    cost_c2 = SCANNERS['C2']['cost']
    cost_c1 = SCANNERS['C1']['cost']
    cost_r1 = SCANNERS['R1']['cost']

    final_coverage_sqm = np.sum(coverage_map)
    final_coverage_ratio = final_coverage_sqm / total_area

    print("\nOptimal Scanner Configuration:")
    print(f"  - C2 Scanners (20m radius): {n_c2}")
    print(f"  - C1 Scanners (5m radius):  {n_c1}")
    print(f"  - R1 Scanners (10m square): {n_r1}")
    
    print(f"\nFinal Achieved Coverage: {final_coverage_sqm} sq.m. ({final_coverage_ratio:.2%})")

    print("\nFinal Cost Calculation:")
    # Here we output each number in the final equation as requested
    print(f"({n_c2} * {cost_c2}) + ({n_c1} * {cost_c1}) + ({n_r1} * {cost_r1}) = {total_cost}")

    print(f"\nOptimal Total Cost: {total_cost}")
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"\n<<<${final_cost}>>>")
