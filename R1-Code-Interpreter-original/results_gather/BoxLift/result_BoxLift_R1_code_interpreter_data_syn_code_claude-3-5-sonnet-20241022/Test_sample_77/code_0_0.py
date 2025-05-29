from itertools import combinations
import copy

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(available_lifters, weight):
    valid_combinations = []
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(range(len(available_lifters)), r):
            lifter_weights = [available_lifters[i] for i in combo]
            if can_lift(lifter_weights, weight):
                valid_combinations.append(list(combo))
    return sorted(valid_combinations, key=len)

def solve_box_lifting():
    boxes = [253, 162, 171, 57, 68, 142, 196, 83, 94, 287, 256, 236, 153, 79, 111, 42]
    lifters = [61, 62, 122, 75, 143, 131]
    
    # Sort boxes from heaviest to lightest
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combos = find_lifter_combinations([lifters[i] for i in available_lifters], box)
            
            if valid_combos:
                # Use the smallest valid combination
                lifter_indices = valid_combos[0]
                actual_lifters = [available_lifters[i] for i in lifter_indices]
                step.append((box, actual_lifters))
                
                # Remove used lifters
                for lifter in sorted(actual_lifters, reverse=True):
                    available_lifters.remove(lifter)
                remaining_boxes.remove(box)
        
        steps.append(step)
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())