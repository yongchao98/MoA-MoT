import itertools

def can_lift(lifters_combination, weight):
    return sum(lifters_combination) >= weight

def find_lifters_for_weight(available_lifters, weight, lifter_capacities):
    # Try single lifter first
    for i, capacity in enumerate(lifter_capacities):
        if i in available_lifters and capacity >= weight:
            return [i]
    
    # Try combinations of lifters
    for size in range(2, len(available_lifters) + 1):
        for combo_indices in itertools.combinations(available_lifters, size):
            combo_capacities = [lifter_capacities[i] for i in combo_indices]
            if sum(combo_capacities) >= weight:
                return list(combo_indices)
    return None

def solve_box_lifting():
    boxes = [87, 123, 302, 308, 370, 86, 182, 79, 135, 170, 198, 243, 45, 212, 183, 174, 335, 134, 58, 131, 124, 360, 192, 245]
    lifters = [167, 156, 74, 177, 77, 60]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            lifters_found = find_lifters_for_weight(available_lifters, box, lifters)
            if lifters_found:
                step.append((box, lifters_found))
                available_lifters -= set(lifters_found)
                remaining_boxes.remove(box)
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(f"<<{result}>>")

solve_box_lifting()