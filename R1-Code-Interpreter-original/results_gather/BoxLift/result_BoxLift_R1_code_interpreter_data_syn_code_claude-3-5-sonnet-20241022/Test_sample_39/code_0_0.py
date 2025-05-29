import itertools

def find_lifter_combinations(box_weight, lifters, used_lifters):
    available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
    
    # Try combinations of different sizes
    for size in range(1, len(available_lifters) + 1):
        for combo in itertools.combinations(available_lifters, size):
            if sum(lifters[i] for i in combo) >= box_weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [256, 53, 286, 205, 261, 209, 242, 143, 354, 55, 72, 161, 268, 305, 43, 246, 146, 331, 84, 128, 88, 358, 71, 93]
    lifters = [189, 58, 58, 79, 123, 192, 116]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign as many boxes as possible in current step
        i = 0
        while i < len(remaining_boxes):
            box_idx, box_weight = remaining_boxes[i]
            lifter_combo = find_lifter_combinations(box_weight, lifters, used_lifters)
            
            if lifter_combo:
                step.append((box_weight, lifter_combo))
                used_lifters.update(lifter_combo)
                remaining_boxes.pop(i)
            else:
                i += 1
                
        if step:
            steps.append(step)
        
        if len(steps) > 8:
            return "No solution found within 8 steps"
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())