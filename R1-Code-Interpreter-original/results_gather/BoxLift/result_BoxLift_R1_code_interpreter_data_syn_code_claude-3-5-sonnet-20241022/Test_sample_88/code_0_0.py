import itertools

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in itertools.combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [88, 313, 344, 387, 392, 370, 99, 389, 340, 391, 131, 147, 334, 251, 157, 144, 290, 167, 141, 79]
    lifter_capacities = [148, 193, 54, 146, 106, 53, 88]
    
    # Sort boxes by weight (descending)
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find lifters that can handle this box
            possible_lifters = []
            for r in range(1, len(available_lifters) + 1):
                for combo in itertools.combinations(available_lifters, r):
                    if sum(lifter_capacities[i] for i in combo) >= box_weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step.append((box_weight, possible_lifters))
                available_lifters -= set(possible_lifters)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()