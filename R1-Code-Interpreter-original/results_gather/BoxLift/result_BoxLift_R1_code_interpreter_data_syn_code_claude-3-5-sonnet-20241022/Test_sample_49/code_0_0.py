import itertools

# Input data
boxes = [351, 64, 363, 78, 283, 68, 211, 200, 288, 295, 319, 271, 305, 219, 310, 365, 390, 56, 396, 332, 130, 136, 221, 386]
lifters = [134, 184, 163, 144, 109, 179, 64]

def can_lift(lifter_combination, weight):
    return sum(lifters[i] for i in lifter_combination) >= weight

def find_lifter_combination(weight, available_lifters, max_lifters=4):
    for n in range(1, max_lifters + 1):
        for combo in itertools.combinations(available_lifters, n):
            if can_lift(combo, weight):
                return list(combo)
    return None

def solve_box_lifting():
    # Sort boxes in descending order
    sorted_boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    steps = []
    remaining_boxes = list(sorted_boxes)
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        remaining_boxes_copy = remaining_boxes.copy()
        for box_idx, box_weight in remaining_boxes_copy:
            if not available_lifters:
                break
                
            lifter_combo = find_lifter_combination(box_weight, available_lifters)
            if lifter_combo:
                step.append((box_weight, lifter_combo))
                available_lifters -= set(lifter_combo)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()