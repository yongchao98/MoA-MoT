import copy

def can_lift(weight, lifters, capacities):
    """Check if combination of lifters can lift the weight"""
    if not lifters:
        return False
    return sum(capacities[i] for i in lifters) >= weight

def find_minimal_lifters(weight, available_lifters, capacities):
    """Find minimal combination of lifters that can lift the weight"""
    if not available_lifters:
        return None
    
    # Try combinations from smallest to largest size
    for size in range(1, len(available_lifters) + 1):
        current_sum = 0
        current_lifters = available_lifters[:size]
        for lifter in current_lifters:
            current_sum += capacities[lifter]
        if current_sum >= weight:
            return current_lifters
    return None

def solve_box_lifting():
    boxes = [112, 109, 151, 97, 375, 276, 159, 378, 284, 275, 120, 60, 169, 398, 125, 268, 110, 92, 166, 349]
    capacities = [102, 116, 199, 175, 148, 161, 161]
    
    # Sort boxes by weight (descending)
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    step = 1
    while remaining_boxes and step <= 6:
        step_solution = []
        available_lifters = list(range(len(capacities)))
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            box_idx, weight = remaining_boxes[i]
            
            # Find minimal combination of lifters for current box
            lifters = find_minimal_lifters(weight, available_lifters, capacities)
            
            if lifters:
                step_solution.append((weight, lifters))
                for lifter in lifters:
                    available_lifters.remove(lifter)
                remaining_boxes.pop(i)
            else:
                i += 1
        
        if step_solution:
            solution.append(step_solution)
        step += 1
    
    # Print solution
    if not remaining_boxes:
        result = ""
        for i, step_sol in enumerate(solution, 1):
            result += f"Step {i}: {step_sol}\n"
        print(result.strip())
    else:
        print("Could not find solution within 6 steps")

solve_box_lifting()