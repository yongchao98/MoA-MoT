from itertools import combinations
import copy

def can_lift(lifters_indices, lifters_capacity, weight):
    return sum(lifters_capacity[i] for i in lifters_indices) >= weight

def find_lifter_combinations(box_weight, available_lifters, lifters_capacity):
    # Try single lifters first
    for i in range(len(available_lifters)):
        if available_lifters[i] and lifters_capacity[i] >= box_weight:
            return [i]
    
    # Try combinations of lifters
    for size in range(2, len(lifters_capacity) + 1):
        for comb in combinations(range(len(lifters_capacity)), size):
            if all(available_lifters[i] for i in comb) and can_lift(comb, lifters_capacity, box_weight):
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [34, 92, 289, 35, 140, 174, 244, 136, 80, 252, 136, 181, 92, 235, 61, 286]
    lifters_capacity = [61, 60, 150, 96, 58, 91]
    
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = [True] * len(lifters_capacity)
        
        # Try to assign as many boxes as possible in current step
        i = 0
        while i < len(remaining_boxes):
            box_idx, box_weight = remaining_boxes[i]
            
            # Find lifters for this box
            lifters = find_lifter_combinations(box_weight, available_lifters, lifters_capacity)
            
            if lifters is not None:
                # Mark lifters as unavailable
                for lifter in lifters:
                    available_lifters[lifter] = False
                step.append((box_weight, lifters))
                remaining_boxes.pop(i)
            else:
                i += 1
        
        if step:
            solution.append(step)
        
        if len(solution) > 6:
            return None
    
    # Format solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

print(solve_box_lifting())