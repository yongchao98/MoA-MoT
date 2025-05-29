import itertools

def find_lifter_combinations(box_weight, lifters, used_lifters):
    available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
    for r in range(1, len(available_lifters) + 1):
        for comb in itertools.combinations(available_lifters, r):
            if sum(lifters[i] for i in comb) >= box_weight:
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [364, 138, 172, 200, 44, 194, 195, 282, 189, 400, 233, 307, 71, 300, 306, 346, 61, 136, 168, 392, 97, 367, 118, 354]
    lifters = [54, 124, 186, 103, 122, 112]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(steps) < 10:
        step = []
        used_lifters = set()
        
        # Try to assign boxes to lifters in current step
        i = 0
        while i < len(remaining_boxes):
            box = remaining_boxes[i]
            lifter_comb = find_lifter_combinations(box, lifters, used_lifters)
            
            if lifter_comb:
                step.append((box, lifter_comb))
                used_lifters.update(lifter_comb)
                remaining_boxes.pop(i)
            else:
                i += 1
                
        if step:
            steps.append(step)
        else:
            break
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    if not remaining_boxes:
        print(result.strip())
    else:
        print("Could not lift all boxes within 10 steps!")

solve_box_lifting()