from collections import defaultdict

def can_lift(box_weight, lifters, used_lifters):
    # Try all possible combinations of available lifters
    n = len(lifters)
    for mask in range(1 << n):
        total_capacity = 0
        lifter_indices = []
        for i in range(n):
            if (mask >> i) & 1:
                if i in used_lifters:
                    break
                total_capacity += lifters[i]
                lifter_indices.append(i)
        else:
            if total_capacity >= box_weight:
                return lifter_indices
    return None

def solve_box_lifting(boxes, lifters):
    boxes = list(enumerate(boxes))  # Keep track of original indices
    steps = []
    
    while boxes:
        step = []
        used_lifters = set()
        
        # Sort remaining boxes by weight (descending)
        boxes.sort(key=lambda x: x[1], reverse=True)
        
        # Try to assign boxes in current step
        remaining_boxes = []
        for box_idx, box_weight in boxes:
            lifter_combination = can_lift(box_weight, lifters, used_lifters)
            if lifter_combination is not None and len(used_lifters) + len(lifter_combination) <= len(lifters):
                step.append((box_weight, lifter_combination))
                used_lifters.update(lifter_combination)
            else:
                remaining_boxes.append((box_idx, box_weight))
        
        if not step:  # If we couldn't assign any boxes in this step
            return None
        
        steps.append(step)
        boxes = remaining_boxes
        
        if len(steps) > 5:  # Check if we exceed maximum steps
            return None
    
    return steps

# Input data
boxes = [71, 55, 11, 10, 71, 91, 14, 69]
lifters = [46, 54, 47]

# Solve and format output
result = solve_box_lifting(boxes, lifters)

if result:
    formatted_output = ""
    for i, step in enumerate(result, 1):
        formatted_output += f"Step {i}: {step}\n"
    print(formatted_output.strip())
else:
    print("No solution found within 5 steps")