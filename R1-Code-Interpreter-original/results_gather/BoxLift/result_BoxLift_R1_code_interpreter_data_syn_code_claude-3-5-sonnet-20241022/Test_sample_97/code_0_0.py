from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Returns list of minimal combinations of lifters that can lift the weight
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                valid_combinations.append(list(combo))
    return sorted(valid_combinations, key=len)

def find_solution():
    boxes = [184, 93, 275, 216, 137, 181, 31, 79, 56, 138, 81, 205, 108, 193, 230, 252]
    lifters = [47, 157, 156, 45, 151, 51]
    
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index[:]
    
    while remaining_boxes and len(solution) < 6:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combos = [c for c in can_lift(lifters, box_weight) 
                          if all(i in available_lifters for i in c)]
            
            if valid_combos:
                best_combo = valid_combos[0]  # Use minimal combination
                step.append((box_weight, best_combo))
                available_lifters -= set(best_combo)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        else:
            return None
    
    return solution if not remaining_boxes else None

# Find and format solution
solution = find_solution()
if solution:
    formatted = []
    for i, step in enumerate(solution, 1):
        step_str = f"Step {i}: {step}"
        formatted.append(step_str)
    print("\n".join(formatted))
else:
    print("No solution found within 6 steps")