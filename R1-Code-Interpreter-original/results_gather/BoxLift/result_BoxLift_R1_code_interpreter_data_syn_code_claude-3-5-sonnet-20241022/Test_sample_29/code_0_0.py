import itertools

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(weight, lifters, used_lifters):
    available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
    
    # Try combinations from smallest to largest size
    for size in range(1, len(available_lifters) + 1):
        for combo in itertools.combinations(available_lifters, size):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting(boxes, lifters):
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign as many boxes as possible in current step
        for box in remaining_boxes[:]:
            lifter_combo = find_lifter_combinations(box, lifters, used_lifters)
            if lifter_combo:
                step.append((box, lifter_combo))
                used_lifters.update(lifter_combo)
                remaining_boxes.remove(box)
        
        if step:
            steps.append(step)
        else:
            break
            
    # Format output
    result = []
    for i, step in enumerate(steps, 1):
        result.append(f"Step {i}: {step}")
    
    return "\n".join(result)

# Input data
boxes = [100, 225, 68, 49, 127, 244, 299, 344, 330, 310, 299, 189, 124, 216, 350, 367, 85, 231, 165, 144, 323, 342, 368, 179]
lifters = [64, 93, 158, 137, 107, 95, 173]

print(solve_box_lifting(boxes, lifters))