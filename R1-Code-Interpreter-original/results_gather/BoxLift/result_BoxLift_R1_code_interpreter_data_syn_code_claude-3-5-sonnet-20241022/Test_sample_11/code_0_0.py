from itertools import combinations
import copy

def can_lift_box(lifters, box_weight):
    # Returns list of lifter combinations that can lift the box
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for comb in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in comb) >= box_weight:
                valid_combinations.append(list(comb))
    return valid_combinations

def solve_box_lifting():
    boxes = [196, 90, 142, 22, 143, 101, 43, 166, 79, 129, 77, 187]
    lifters = [71, 58, 104, 105]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    boxes = [(i, w) for i, w in boxes]  # Keep original index
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 6:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combs = [c for c in can_lift_box([lifters[i] for i in available_lifters], box_weight)
                         if all(i in available_lifters for i in c)]
            
            if valid_combs:
                # Use the combination with minimum number of lifters
                best_comb = min(valid_combs, key=len)
                step.append((box_weight, best_comb))
                available_lifters -= set(best_comb)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        else:
            break
    
    # Format and print solution
    if not remaining_boxes and len(solution) <= 6:
        result = []
        for i, step in enumerate(solution, 1):
            step_str = f"Step {i}: {step}"
            result.append(step_str)
        print("<<<" + "\n".join(result) + ">>>")
    else:
        print("<<<No solution found within 6 steps>>>")

solve_box_lifting()