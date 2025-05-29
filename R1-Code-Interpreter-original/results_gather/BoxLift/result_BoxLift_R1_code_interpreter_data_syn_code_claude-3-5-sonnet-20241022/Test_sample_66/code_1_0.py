from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [181, 52, 82, 130, 92, 195, 101, 106, 189, 42, 196, 46]
    lifter_capacities = [45, 30, 52, 104, 30]
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifter_capacities)))
        
        # Sort boxes in descending order for greedy approach
        remaining_boxes.sort(reverse=True)
        boxes_handled_this_step = []
        
        # Try to assign boxes to lifters
        for i, box in enumerate(remaining_boxes):
            if not available_lifters:  # No more lifters available for this step
                break
                
            lifter_caps = [lifter_capacities[j] for j in available_lifters]
            lifters_needed = can_lift(lifter_caps, box)
            
            if lifters_needed is not None:
                # Convert relative indices to actual available lifter indices
                actual_lifters = [available_lifters[j] for j in lifters_needed]
                step.append((box, actual_lifters))
                boxes_handled_this_step.append(box)
                
                # Remove used lifters (in reverse order to maintain correct indices)
                for idx in sorted(lifters_needed, reverse=True):
                    available_lifters.pop(idx)
        
        # Remove handled boxes
        for box in boxes_handled_this_step:
            remaining_boxes.remove(box)
            
        if step:
            steps.append(step)
        else:
            break
            
        if len(steps) > 8:
            break
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    if not remaining_boxes and len(steps) <= 8:
        print(result.strip())
    else:
        print("No solution found within 8 steps")

solve_box_lifting()