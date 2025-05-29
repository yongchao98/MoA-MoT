from itertools import combinations
import copy

def can_lift_box(lifters_capacity, box_weight):
    if sum(lifters_capacity) < box_weight:
        return False
    if max(lifters_capacity) >= box_weight:
        return True
    
    for r in range(2, len(lifters_capacity) + 1):
        for comb in combinations(range(len(lifters_capacity)), r):
            if sum(lifters_capacity[i] for i in comb) >= box_weight:
                return True
    return False

def find_lifters_for_box(available_lifters, lifters_capacity, box_weight):
    # Try single lifter first
    for i in available_lifters:
        if lifters_capacity[i] >= box_weight:
            return [i]
    
    # Try combinations of lifters
    for r in range(2, len(available_lifters) + 1):
        for comb in combinations(available_lifters, r):
            if sum(lifters_capacity[i] for i in comb) >= box_weight:
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
    lifters = [164, 182, 169, 127, 118, 60]
    
    # Sort boxes from heaviest to lightest
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in this step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            lifters_for_box = find_lifters_for_box(list(available_lifters), lifters, box_weight)
            if lifters_for_box:
                step.append((box_weight, lifters_for_box))
                available_lifters -= set(lifters_for_box)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
        else:
            break
            
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(f"<<{result}>>")

solve_box_lifting()