from itertools import combinations
import copy

def can_lift(lifters_combo, weight):
    return sum(lifters_combo) >= weight

def find_lifter_combinations(box, available_lifters, lifter_capacities):
    # Try combinations of lifters from smallest to largest size
    lifter_indices = list(range(len(available_lifters)))
    for size in range(1, len(available_lifters) + 1):
        for combo_indices in combinations(lifter_indices, size):
            # Check if these lifters are all available and can lift the box
            if all(available_lifters[i] for i in combo_indices):
                combo_capacity = sum(lifter_capacities[i] for i in combo_indices)
                if combo_capacity >= box:
                    return list(combo_indices)
    return None

def solve_box_lifting():
    boxes = [257, 47, 186, 194, 237, 135, 72, 204, 268, 83, 242, 41, 214, 157, 115, 242]
    lifter_capacities = [77, 41, 40, 54, 132]
    
    # Sort boxes in descending order
    boxes_with_indices = list(enumerate(boxes))
    boxes_with_indices.sort(key=lambda x: x[1], reverse=True)
    
    steps = []
    remaining_boxes = boxes_with_indices.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = [True] * len(lifter_capacities)
        
        # Try to assign as many boxes as possible in this step
        for box_idx, box_weight in remaining_boxes[:]:
            lifters = find_lifter_combinations(box_weight, available_lifters, lifter_capacities)
            if lifters:
                # Mark these lifters as unavailable for this step
                for lifter in lifters:
                    available_lifters[lifter] = False
                step.append((box_weight, lifters))
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
        else:
            break
            
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()