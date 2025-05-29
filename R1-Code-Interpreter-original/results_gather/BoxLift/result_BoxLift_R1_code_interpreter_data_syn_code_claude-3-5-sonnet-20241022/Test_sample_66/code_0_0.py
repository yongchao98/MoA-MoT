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
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            box = remaining_boxes[i]
            lifters_needed = can_lift([lifter_capacities[j] for j in available_lifters], box)
            
            if lifters_needed is not None:
                # Convert lifters_needed indices to actual available lifter indices
                actual_lifters = [available_lifters[j] for j in lifters_needed]
                step.append((box, actual_lifters))
                # Remove used lifters
                for lifter in lifters_needed:
                    available_lifters.pop(lifter)
                remaining_boxes.pop(i)
            else:
                i += 1
        
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