from itertools import combinations
import copy

def can_lift_box(lifters_capacity, box_weight):
    if not lifters_capacity:
        return None
    # Try combinations of lifters from 1 to len(lifters)
    for i in range(1, len(lifters_capacity) + 1):
        for comb in combinations(range(len(lifters_capacity)), i):
            if sum(lifters_capacity[j] for j in comb) >= box_weight:
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [280, 109, 195, 129, 54, 216, 157, 193, 44, 233, 171, 218, 204, 120, 131, 151]
    lifters = [49, 131, 66, 124, 91]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifters)))
        
        # Try to assign boxes to lifters
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            lifters_capacity = [lifters[i] for i in available_lifters]
            lifter_combination = can_lift_box(lifters_capacity, box_weight)
            
            if lifter_combination:
                # Convert relative indices to absolute indices
                actual_lifters = [available_lifters[i] for i in lifter_combination]
                step.append((box_weight, actual_lifters))
                remaining_boxes.remove((box_idx, box_weight))
                # Remove used lifters
                for idx in sorted(lifter_combination, reverse=True):
                    available_lifters.pop(idx)
        
        if step:
            steps.append(step)
        
        if len(steps) > 8:
            return None
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

print(solve_box_lifting())