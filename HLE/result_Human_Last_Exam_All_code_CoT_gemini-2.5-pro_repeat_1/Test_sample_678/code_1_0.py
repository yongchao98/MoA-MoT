import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Problem Definition
    ROOM_W, ROOM_H = 140, 110
    TARGET_COVERAGE_RATIO = 0.88
    PLACEMENT_STEP = 5

    SCANNERS = {
        'C2': {'shape': 'circle', 'param': 20, 'cost': 20000}, # param is radius
        'C1': {'shape': 'circle', 'param': 5,  'cost': 1600},  # param is radius (from 10m diameter)
        'R1': {'shape': 'square', 'param': 10, 'cost': 2000}   # param is side length
    }

    # 2. Discretization Setup
    # Create a grid of points representing the room area (1 point per m^2)
    X_grid, Y_grid = np.mgrid[0:ROOM_W + 1, 0:ROOM_H + 1]
    total_points = (ROOM_W + 1) * (ROOM_H + 1)
    target_covered_points = int(total_points * TARGET_COVERAGE_RATIO)

    # Generate all possible center coordinates for the scanners
    possible_centers = []
    for x in range(0, ROOM_W + 1, PLACEMENT_STEP):
        for y in range(0, ROOM_H + 1, PLACEMENT_STEP):
            possible_centers.append((x, y))

    # Helper function to calculate the coverage mask for a given scanner placement
    def get_coverage_mask(scanner_type_info, center):
        cx, cy = center
        shape = scanner_type_info['shape']
        param = scanner_type_info['param']
        
        if shape == 'circle':
            radius = param
            return ((X_grid - cx)**2 + (Y_grid - cy)**2) <= radius**2
        elif shape == 'square':
            half_side = param / 2
            return (X_grid >= cx - half_side) & (X_grid <= cx + half_side) & \
                   (Y_grid >= cy - half_side) & (Y_grid <= cy + half_side)
        return np.zeros_like(X_grid, dtype=bool)

    # 3. Greedy Algorithm Implementation
    uncovered_mask = np.ones_like(X_grid, dtype=bool)
    total_cost = 0
    num_covered_points = 0
    selected_scanners_summary = {}

    print(f"Room Area: {ROOM_W}x{ROOM_H} m^2. Target Coverage: {TARGET_COVERAGE_RATIO:.0%}")
    print(f"Discretized into {total_points} points. Need to cover {target_covered_points} points.")
    print("Starting optimization...\n")

    iteration = 0
    while num_covered_points < target_covered_points:
        iteration += 1
        best_candidate = None
        max_effectiveness_score = -1
        best_newly_covered_mask = None
        
        # Find the most cost-effective scanner placement in this iteration
        for scanner_name, scanner_info in SCANNERS.items():
            cost = scanner_info['cost']
            for center in possible_centers:
                coverage_mask = get_coverage_mask(scanner_info, center)
                newly_covered_mask = coverage_mask & uncovered_mask
                num_newly_covered = np.sum(newly_covered_mask)
                
                if num_newly_covered > 0:
                    effectiveness_score = num_newly_covered / cost
                    if effectiveness_score > max_effectiveness_score:
                        max_effectiveness_score = effectiveness_score
                        best_candidate = {
                            'name': scanner_name,
                            'cost': cost,
                            'newly_covered': num_newly_covered
                        }
                        best_newly_covered_mask = newly_covered_mask

        if best_candidate is None:
            print("No more points can be covered. Stopping.")
            break
            
        # Add the best scanner from this iteration to our solution
        total_cost += best_candidate['cost']
        num_covered_points += best_candidate['newly_covered']
        uncovered_mask[best_newly_covered_mask] = False
        
        # Update summary
        name = best_candidate['name']
        selected_scanners_summary[name] = selected_scanners_summary.get(name, 0) + 1
        
        print(f"Step {iteration}: Added {best_candidate['name']} scanner. "
              f"Coverage: {num_covered_points / total_points:.2%}. "
              f"Total Cost: ${total_cost}")

    # 4. Output Final Results
    print("\n--- Optimization Finished ---")
    print(f"Final Coverage: {num_covered_points / total_points:.2%} ({num_covered_points}/{total_points} points)")
    
    print("\n--- Final Cost Calculation ---")
    equation_parts = []
    # Ensure a consistent order for the equation (C2, R1, C1)
    for name in ['C2', 'R1', 'C1']:
        if name in selected_scanners_summary:
            count = selected_scanners_summary[name]
            cost = SCANNERS[name]['cost']
            equation_parts.append(f"{count} * {cost}")

    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {total_cost}")

    print("\nThe optimal total cost is:")
    print(f"<<<{total_cost}>>>")

# Run the solver
solve_scanner_placement()