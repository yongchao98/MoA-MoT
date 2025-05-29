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
    boxes = [45, 32, 84, 101, 161, 118, 30, 152, 159, 20, 110, 103]
    lifter_capacities = [39, 32, 112, 102, 75]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(steps) < 5:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters
        current_boxes = remaining_boxes.copy()
        for box in current_boxes:
            if not available_lifters:
                break
                
            lifters = find_lifter_combinations(box, available_lifters, lifter_capacities)
            if lifters:
                step.append((box, lifters))
                remaining_boxes.remove(box)
                available_lifters -= set(lifters)
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()