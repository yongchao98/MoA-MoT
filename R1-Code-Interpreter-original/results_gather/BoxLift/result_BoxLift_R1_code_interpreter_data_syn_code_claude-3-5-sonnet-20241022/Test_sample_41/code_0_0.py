import itertools

def can_lift(box_weight, lifters_combination, lifter_capacities):
    return sum(lifter_capacities[i] for i in lifters_combination) >= box_weight

def find_solution():
    boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
    lifter_capacities = [95, 111, 115, 35]
    n_lifters = len(lifter_capacities)
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes and len(steps) < 5:
        step = []
        available_lifters = set(range(n_lifters))
        
        # Try to assign boxes to lifters
        remaining_boxes_copy = remaining_boxes.copy()
        for box in remaining_boxes_copy:
            if not available_lifters:
                break
                
            # Try combinations of available lifters
            assigned = False
            for size in range(1, len(available_lifters) + 1):
                for lifters_combo in itertools.combinations(available_lifters, size):
                    if can_lift(box, lifters_combo, lifter_capacities):
                        step.append((box, list(lifters_combo)))
                        available_lifters -= set(lifters_combo)
                        remaining_boxes.remove(box)
                        assigned = True
                        break
                if assigned:
                    break
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(f"<<<{result.strip()}>>>")

find_solution()