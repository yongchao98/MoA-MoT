import copy

def can_lift_box(box_weight, lifters, used_lifters):
    # Try single lifter first
    for i, capacity in enumerate(lifters):
        if i not in used_lifters and capacity >= box_weight:
            return [i]
    
    # Try combinations of lifters
    total_lifters = len(lifters)
    best_combination = None
    min_lifters_needed = float('inf')
    
    # Try different combinations of lifters
    for i in range(total_lifters):
        if i in used_lifters:
            continue
        current_sum = lifters[i]
        current_lifters = [i]
        
        for j in range(i + 1, total_lifters):
            if j in used_lifters:
                continue
            if current_sum + lifters[j] >= box_weight:
                current_lifters.append(j)
                if len(current_lifters) < min_lifters_needed:
                    min_lifters_needed = len(current_lifters)
                    best_combination = current_lifters.copy()
                break
            current_sum += lifters[j]
            current_lifters.append(j)
    
    return best_combination

def solve_box_lifting():
    boxes = [166, 144, 53, 213, 51, 156, 197, 311, 177, 358, 172, 134, 179, 145, 91, 188, 352, 294, 292, 88, 97, 394, 123, 294]
    lifters = [88, 185, 145, 195, 147, 145]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign as many boxes as possible in current step
        remaining_boxes_copy = remaining_boxes.copy()
        for box in remaining_boxes_copy:
            lifter_combination = can_lift_box(box, lifters, used_lifters)
            if lifter_combination:
                step.append((box, lifter_combination))
                used_lifters.update(lifter_combination)
                remaining_boxes.remove(box)
            
            if len(used_lifters) == len(lifters):
                break
        
        steps.append(step)
    
    # Print result
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(f"<<{result}>>")

solve_box_lifting()