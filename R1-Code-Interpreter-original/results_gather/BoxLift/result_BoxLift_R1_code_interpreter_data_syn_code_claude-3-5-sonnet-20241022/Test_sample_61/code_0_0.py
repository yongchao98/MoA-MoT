from itertools import combinations
import copy

def can_lift_box(lifters, box_weight):
    # Returns list of lifters that can lift the box (single or combined)
    all_combinations = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= box_weight:
                all_combinations.append(list(combo))
    return all_combinations

def find_solution():
    boxes = [128, 144, 42, 49, 132, 118, 32, 48, 21, 94, 85, 51]
    lifter_capacities = [85, 90, 40, 102]
    
    # Sort boxes from heaviest to lightest
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    boxes = [(i, w) for i, w in boxes]  # Keep original index
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 5:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find all possible lifter combinations for this box
            possible_lifters = [combo for combo in can_lift_box(
                [lifter_capacities[i] for i in available_lifters], 
                box_weight
            ) if all(l in available_lifters for l in combo)]
            
            if possible_lifters:
                # Use the combination with minimum number of lifters
                best_combo = min(possible_lifters, key=len)
                step.append((box_weight, best_combo))
                for lifter in best_combo:
                    available_lifters.remove(lifter)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        else:
            break
    
    # Format the solution
    if not remaining_boxes and len(solution) <= 5:
        result = []
        for i, step in enumerate(solution, 1):
            result.append(f"Step {i}: {step}")
        print("<<<" + "\n".join(result) + ">>>")
    else:
        print("<<<No solution found within 5 steps>>>")

find_solution()