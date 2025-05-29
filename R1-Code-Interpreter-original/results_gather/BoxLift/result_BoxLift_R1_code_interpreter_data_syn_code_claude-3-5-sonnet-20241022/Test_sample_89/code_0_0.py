from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [227, 106, 388, 173, 63, 178, 261, 205, 323, 124, 386, 81, 71, 127, 134, 212, 150, 114, 41, 277, 123, 152, 47, 47]
    lifter_capacities = [166, 106, 122, 181, 102, 119, 147]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes in current step
        current_boxes = remaining_boxes.copy()
        for box in current_boxes:
            if not available_lifters:
                break
                
            # Find lifters that can handle this box
            lifters_for_box = []
            test_lifters = [lifter_capacities[i] for i in available_lifters]
            test_lifter_indices = list(available_lifters)
            
            lifter_combo = can_lift([lifter_capacities[i] for i in available_lifters], box)
            
            if lifter_combo:
                # Convert relative indices to absolute indices
                actual_lifters = [test_lifter_indices[i] for i in lifter_combo]
                step.append((box, actual_lifters))
                remaining_boxes.remove(box)
                available_lifters -= set(actual_lifters)
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())