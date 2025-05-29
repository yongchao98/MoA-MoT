import copy

def can_lift_box(lifters_available, lifter_capacities, box_weight):
    # Find combination of available lifters that can lift the box
    if not lifters_available:
        return None
    
    # Try single lifter first
    for i in lifters_available:
        if lifter_capacities[i] >= box_weight:
            return [i]
    
    # Try combinations of lifters
    total_lifters = len(lifters_available)
    for size in range(2, len(lifters_available) + 1):
        for i in range(total_lifters):
            for j in range(i + 1, total_lifters):
                lifters = [lifters_available[i], lifters_available[j]]
                if sum(lifter_capacities[k] for k in lifters) >= box_weight:
                    return lifters
    return None

def solve_box_lifting():
    boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
    lifter_capacities = [122, 116, 51, 165, 95, 150, 179]
    
    # Sort boxes by weight (descending)
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        lifters_available = list(range(len(lifter_capacities)))
        
        # Try to assign as many boxes as possible in this step
        i = 0
        while i < len(remaining_boxes):
            box_idx, box_weight = remaining_boxes[i]
            lifters = can_lift_box(lifters_available, lifter_capacities, box_weight)
            
            if lifters:
                step.append((box_weight, lifters))
                # Remove used lifters
                for lifter in lifters:
                    lifters_available.remove(lifter)
                remaining_boxes.pop(i)
            else:
                i += 1
                
        solution.append(step)
        
        if len(solution) > 7:
            print("Cannot solve within 7 steps!")
            return None
    
    # Format and print solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    print(result.strip())

solve_box_lifting()