import itertools

def can_lift(weight, lifters, capacities):
    return sum(capacities[i] for i in lifters) >= weight

def find_solution():
    boxes = [64, 70, 386, 351, 113, 77, 314, 333, 266, 399, 193, 44, 181, 200, 238, 175, 370, 118, 337, 134]
    capacities = [140, 115, 159, 147, 129, 112]
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)  # Sort boxes by weight
    n_lifters = len(capacities)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(n_lifters))
        
        # Try to assign heaviest boxes first
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            for n_lifters_needed in range(1, len(available_lifters) + 1):
                found = False
                for lifter_combo in itertools.combinations(available_lifters, n_lifters_needed):
                    if can_lift(box_weight, lifter_combo, capacities):
                        step.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.remove((box_idx, box_weight))
                        found = True
                        break
                if found:
                    break
        
        if step:
            steps.append(step)
        if len(steps) > 7:
            return None
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

solution = find_solution()
print(solution)