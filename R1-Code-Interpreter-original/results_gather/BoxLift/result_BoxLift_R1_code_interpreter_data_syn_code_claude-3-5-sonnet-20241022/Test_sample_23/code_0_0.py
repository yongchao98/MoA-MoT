import itertools

def can_lift(lifters, weight):
    # Try all possible combinations of lifters to find if they can lift the weight
    for r in range(1, len(lifters) + 1):
        for combo in itertools.combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [254, 205, 272, 360, 320, 328, 267, 223, 337, 238, 134, 235, 388, 160, 244, 378, 349, 262, 104, 396]
    lifters = [132, 86, 81, 199, 93, 92]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box_weight, box_idx in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find minimum combination of available lifters that can lift this box
            possible_lifters = [i for i in available_lifters]
            lifter_combo = can_lift([lifters[i] for i in possible_lifters], box_weight)
            
            if lifter_combo:
                # Convert local lifter indices to global indices
                global_lifter_combo = [possible_lifters[i] for i in lifter_combo]
                step.append((box_weight, global_lifter_combo))
                remaining_boxes.remove((box_weight, box_idx))
                available_lifters -= set(global_lifter_combo)
        
        if step:
            steps.append(step)
    
    # Format and print the solution
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    print(result.strip())

solve_box_lifting()