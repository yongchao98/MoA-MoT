from itertools import combinations

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifters_for_weight(available_lifters, weight, lifter_indices):
    # Try combinations of lifters from smallest to largest groups
    for i in range(1, len(available_lifters) + 1):
        for combo in combinations(range(len(available_lifters)), i):
            lifter_weights = [available_lifters[j] for j in combo]
            if can_lift(lifter_weights, weight):
                return [lifter_indices[j] for j in combo]
    return None

def solve_box_lifting():
    boxes = [202, 297, 200, 241, 156, 177, 70, 78, 133, 176, 75, 70, 167, 240, 299, 91]
    lifter_capacities = [56, 160, 118, 138, 160, 150]
    
    # Sort boxes from heaviest to lightest
    boxes = sorted(boxes, reverse=True)
    steps = []
    
    while boxes:
        step = []
        available_lifters = lifter_capacities.copy()
        lifter_indices = list(range(len(lifter_capacities)))
        
        # Try to assign as many boxes as possible in this step
        i = 0
        while i < len(boxes):
            box = boxes[i]
            lifters = find_lifters_for_weight(available_lifters, box, lifter_indices)
            
            if lifters:
                # Remove used lifters
                for lifter_idx in lifters:
                    idx = lifter_indices.index(lifter_idx)
                    available_lifters.pop(idx)
                    lifter_indices.pop(idx)
                
                step.append((box, lifters))
                boxes.pop(i)
            else:
                i += 1
                
        if step:
            steps.append(step)
            
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()