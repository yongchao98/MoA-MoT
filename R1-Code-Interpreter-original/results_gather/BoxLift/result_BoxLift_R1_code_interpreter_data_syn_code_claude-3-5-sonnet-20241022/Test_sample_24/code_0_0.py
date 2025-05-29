import copy

def can_lift(weight, lifters_capacity, used_lifters):
    # Try different combinations of available lifters
    n = len(lifters_capacity)
    available = [i for i in range(n) if i not in used_lifters]
    
    # Try single lifter first
    for i in available:
        if lifters_capacity[i] >= weight:
            return [i]
    
    # Try combinations of lifters
    for size in range(2, len(available) + 1):
        for i in range(len(available) - size + 1):
            lifters = available[i:i+size]
            if sum(lifters_capacity[j] for j in lifters) >= weight:
                return lifters
    return None

def solve_box_lifting():
    boxes = [119, 110, 359, 148, 237, 71, 309, 213, 117, 394, 79, 353, 230, 77, 242, 239, 66, 259, 212, 358, 194, 397, 274, 101]
    lifters = [60, 103, 156, 108, 179, 110]
    
    # Sort boxes by weight (descending)
    boxes = sorted(boxes, reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign as many boxes as possible in current step
        remaining_boxes_copy = remaining_boxes.copy()
        for box in remaining_boxes_copy:
            lifter_combination = can_lift(box, lifters, used_lifters)
            if lifter_combination:
                step.append((box, lifter_combination))
                used_lifters.update(lifter_combination)
                remaining_boxes.remove(box)
            
            if len(used_lifters) == len(lifters):
                break
        
        solution.append(step)
    
    # Format and print solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()