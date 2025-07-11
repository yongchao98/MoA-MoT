import math
from collections import defaultdict

def solve_scanner_placement():
    """
    Formulates the museum scanner problem as a set cover problem and solves it
    using a greedy algorithm to find a near-optimal solution.
    """
    # 1. Define problem constants
    ROOM_W, ROOM_H = 140, 110
    PLACEMENT_STEP = 5
    COVERAGE_TARGET_PERCENT = 0.88

    TOTAL_AREA = ROOM_W * ROOM_H
    TARGET_COVERAGE_AREA = math.floor(TOTAL_AREA * COVERAGE_TARGET_PERCENT)

    SCANNERS = {
        'C2': {'shape': 'circle', 'dim': 20, 'cost': 20000}, # dim is radius
        'C1': {'shape': 'circle', 'dim': 5,  'cost': 1600},  # dim is radius (10m diameter)
        'R1': {'shape': 'square', 'dim': 10, 'cost': 2000}   # dim is side
    }

    print("Problem Setup:")
    print(f"Room Dimensions: {ROOM_W}m x {ROOM_H}m (Total Area: {TOTAL_AREA} m^2)")
    print(f"Required Coverage: {COVERAGE_TARGET_PERCENT*100}% ({TARGET_COVERAGE_AREA} m^2)\n")

    # 2. Pre-calculate footprints for all possible scanner placements
    print("Pre-calculating scanner coverage footprints...")
    potential_scanners = []
    center_x_coords = range(0, ROOM_W + 1, PLACEMENT_STEP)
    center_y_coords = range(0, ROOM_H + 1, PLACEMENT_STEP)

    for name, props in SCANNERS.items():
        for cx in center_x_coords:
            for cy in center_y_coords:
                footprint = set()
                if props['shape'] == 'circle':
                    radius = props['dim']
                    # Iterate only over the bounding box of the circle for efficiency
                    for x in range(max(0, cx - radius), min(ROOM_W, cx + radius + 1)):
                        for y in range(max(0, cy - radius), min(ROOM_H, cy + radius + 1)):
                            if (x - cx)**2 + (y - cy)**2 <= radius**2:
                                footprint.add((x, y))
                elif props['shape'] == 'square':
                    side = props['dim']
                    half_side = side / 2
                    # Iterate only over the bounding box of the square
                    for x in range(max(0, math.ceil(cx - half_side)), min(ROOM_W, math.floor(cx + half_side) + 1)):
                         for y in range(max(0, math.ceil(cy - half_side)), min(ROOM_H, math.floor(cy + half_side) + 1)):
                            if abs(x - cx) <= half_side and abs(y - cy) <= half_side:
                                footprint.add((x, y))
                
                if footprint:
                    potential_scanners.append({
                        'name': name,
                        'cost': props['cost'],
                        'center': (cx, cy),
                        'footprint': footprint
                    })

    print(f"Generated {len(potential_scanners)} potential scanner placements.\n")

    # 3. Greedy algorithm to find the solution
    print("Running greedy algorithm to find the best combination...")
    total_cost = 0
    covered_points = set()
    solution_scanners = []

    while len(covered_points) < TARGET_COVERAGE_AREA:
        best_scanner = None
        max_effectiveness = -1
        best_new_points = set()

        for scanner in potential_scanners:
            # Calculate how many NEW points this scanner covers
            newly_covered = scanner['footprint'] - covered_points
            
            if not newly_covered:
                continue

            # Effectiveness = new area covered / cost
            effectiveness = len(newly_covered) / scanner['cost']
            
            if effectiveness > max_effectiveness:
                max_effectiveness = effectiveness
                best_scanner = scanner
                best_new_points = newly_covered

        if best_scanner is None:
            print("Warning: No more area can be covered.")
            break

        # Add the best scanner found in this iteration to the solution
        solution_scanners.append(best_scanner)
        covered_points.update(best_new_points)
        total_cost += best_scanner['cost']
        
        print(f"  -> Added {best_scanner['name']} at {best_scanner['center']}. "
              f"Coverage: {len(covered_points)}/{TARGET_COVERAGE_AREA} m^2. "
              f"Cost: ${total_cost}")

    print("\n--- Optimization Complete ---")

    # 4. Summarize and print results
    scanner_counts = defaultdict(int)
    for scanner in solution_scanners:
        scanner_counts[scanner['name']] += 1

    cost_c2 = SCANNERS['C2']['cost']
    cost_c1 = SCANNERS['C1']['cost']
    cost_r1 = SCANNERS['R1']['cost']
    
    num_c2 = scanner_counts.get('C2', 0)
    num_c1 = scanner_counts.get('C1', 0)
    num_r1 = scanner_counts.get('R1', 0)

    final_coverage_percent = (len(covered_points) / TOTAL_AREA) * 100

    print("\nFinal Scanner Configuration:")
    print(f"- C2 Scanners (r=20m, ${cost_c2}): {num_c2}")
    print(f"- C1 Scanners (d=10m, ${cost_c1}): {num_c1}")
    print(f"- R1 Scanners (s=10m, ${cost_r1}): {num_r1}")

    print("\nFinal Cost Calculation:")
    cost_str_parts = []
    if num_c2 > 0: cost_str_parts.append(f"{num_c2} * {cost_c2}")
    if num_c1 > 0: cost_str_parts.append(f"{num_c1} * {cost_c1}")
    if num_r1 > 0: cost_str_parts.append(f"{num_r1} * {cost_r1}")
    
    print(f"{' + '.join(cost_str_parts)} = ${total_cost}")

    print(f"\nTotal Achieved Coverage: {len(covered_points)} m^2 ({final_coverage_percent:.2f}%)")
    print(f"Optimal Total Cost: ${total_cost}")
    
    # Final answer block
    print(f"\n<<<{total_cost}>>>")


if __name__ == '__main__':
    solve_scanner_placement()