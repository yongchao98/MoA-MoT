import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Configuration ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    GRID_STEP = 5 # meters for scanner center placement
    TARGET_COVERAGE_FRAC = 0.88

    TOTAL_AREA = ROOM_W * ROOM_H
    TARGET_COVERAGE_AREA = TOTAL_AREA * TARGET_COVERAGE_FRAC

    SCANNERS = [
        {
            "name": "C2",
            "type": "circle",
            "param": 20,  # radius
            "cost": 20000,
        },
        {
            "name": "C1",
            "type": "circle",
            "param": 5,   # radius (from 10m diameter)
            "cost": 1600,
        },
        {
            "name": "R1",
            "type": "square",
            "param": 10,  # side
            "cost": 2000,
        },
    ]

    # --- Initialization ---
    print("Starting optimization...")
    print(f"Room Area: {TOTAL_AREA} m^2")
    print(f"Target Coverage: {TARGET_COVERAGE_AREA:.2f} m^2 ({TARGET_COVERAGE_FRAC:.0%})")
    print("-" * 30)

    # 1. Generate possible scanner locations
    possible_locations = []
    for x in range(0, ROOM_W + 1, GRID_STEP):
        for y in range(0, ROOM_H + 1, GRID_STEP):
            possible_locations.append((x, y))

    # 2. Initialize coverage map (1 point per square meter)
    # 0 = uncovered, 1 = covered.
    coverage_map = np.zeros((ROOM_H, ROOM_W), dtype=np.int8)
    
    # Create coordinate grids once for vectorized calculations
    yy, xx = np.mgrid[:ROOM_H, :ROOM_W]

    # --- Greedy Algorithm ---
    scanners_placed = []
    total_cost = 0
    current_covered_area = 0

    while current_covered_area < TARGET_COVERAGE_AREA:
        best_effectiveness = -1
        best_choice = None

        # Find the most cost-effective next scanner to place
        for scanner in SCANNERS:
            for loc in possible_locations:
                cx, cy = loc
                
                # Generate a mask for the current candidate scanner
                if scanner['type'] == 'circle':
                    r = scanner['param']
                    # Use squared distance to avoid sqrt
                    mask = ((xx - cx)**2 + (yy - cy)**2) <= r**2
                elif scanner['type'] == 'square':
                    s_half = scanner['param'] / 2
                    mask = (np.abs(xx - cx) < s_half) & (np.abs(yy - cy) < s_half)

                # Calculate new coverage (points in mask but not yet covered)
                newly_covered_mask = mask & (coverage_map == 0)
                new_area = np.sum(newly_covered_mask)

                if new_area == 0:
                    continue

                effectiveness = new_area / scanner['cost']

                if effectiveness > best_effectiveness:
                    best_effectiveness = effectiveness
                    best_choice = {
                        "scanner": scanner,
                        "location": loc,
                        "new_area": new_area,
                        "mask": newly_covered_mask
                    }

        if best_choice is None:
            print("No further improvement possible. Stopping.")
            break

        # Place the best scanner found in this iteration
        chosen_scanner = best_choice['scanner']
        chosen_loc = best_choice['location']
        new_area_gain = best_choice['new_area']
        
        scanners_placed.append({
            "name": chosen_scanner['name'],
            "location": chosen_loc,
            "cost": chosen_scanner['cost']
        })
        total_cost += chosen_scanner['cost']
        coverage_map[best_choice['mask']] = 1 # Update coverage map
        current_covered_area = np.sum(coverage_map)

        print(f"Placed {chosen_scanner['name']} at {chosen_loc}. "
              f"Cost: {chosen_scanner['cost']}. "
              f"New Area: {new_area_gain}. "
              f"Total Cost: {total_cost}. "
              f"Coverage: {current_covered_area / TOTAL_AREA:.2%}")

    # --- Final Output ---
    print("\n--- Optimization Finished ---")
    print(f"Final Cost: {total_cost}")
    print(f"Final Coverage: {current_covered_area / TOTAL_AREA:.4f} ({current_covered_area} / {TOTAL_AREA} m^2)")
    print(f"Number of scanners: {len(scanners_placed)}")

    # Summarize scanner counts and build the final equation
    counts = {}
    for s in scanners_placed:
        counts[s['name']] = counts.get(s['name'], 0) + 1
    
    print("\nScanners used:")
    cost_breakdown_parts = []
    for name, count in sorted(counts.items()):
        cost_per_scanner = next(s['cost'] for s in SCANNERS if s['name'] == name)
        print(f"- {name}: {count}")
        cost_breakdown_parts.append(f"{count} * {cost_per_scanner}")

    final_equation = " + ".join(cost_breakdown_parts)
    print("\nFinal Cost Calculation:")
    print(f"{final_equation} = {total_cost}")
    
    print(f"<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()