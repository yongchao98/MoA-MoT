from itertools import combinations
import copy

def can_lift(weight, lifters_capacity):
    return sum(lifters_capacity) >= weight

def find_lifter_combination(weight, available_lifters, lifter_capacities):
    for r in range(1, len(available_lifters) + 1):
        for comb in combinations(available_lifters, r):
            total_capacity = sum(lifter_capacities[i] for i in comb)
            if total_capacity >= weight:
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [247, 236, 255, 64, 58, 144, 218, 200, 320, 87, 256, 111, 379, 318, 367, 88, 136, 125, 172, 93]
    lifter_capacities = [63, 51, 130, 110, 117, 132]
    
    # Sort boxes by weight (descending)
    boxes = sorted(boxes, reverse=True)
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifter_capacities)))
        
        # Try to assign boxes in current step
        boxes_to_remove = []
        for box in remaining_boxes:
            if not available_lifters:
                break
                
            lifters = find_lifter_combination(box, available_lifters, lifter_capacities)
            if lifters:
                step.append((box, lifters))
                for lifter in lifters:
                    available_lifters.remove(lifter)
                boxes_to_remove.append(box)
        
        # Remove assigned boxes
        for box in boxes_to_remove:
            remaining_boxes.remove(box)
            
        steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()