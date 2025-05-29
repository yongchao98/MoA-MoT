import itertools

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(box_weight, available_lifters, lifter_capacities):
    # Try combinations from smallest to largest size
    for i in range(1, len(available_lifters) + 1):
        for combo in itertools.combinations(available_lifters, i):
            if sum(lifter_capacities[j] for j in combo) >= box_weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [49, 239, 263, 284, 45, 159, 306, 398, 279, 137, 218, 387, 161, 287, 171, 345, 373, 271, 60, 396]
    lifter_capacities = [151, 91, 70, 108, 182, 200]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifter_capacities)))
        
        # Try to assign as many boxes as possible in this step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            lifters = find_lifter_combinations(box, available_lifters, lifter_capacities)
            if lifters:
                step.append((box, lifters))
                remaining_boxes.remove(box)
                for lifter in lifters:
                    available_lifters.remove(lifter)
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()