import itertools

def can_lift_box(box_weight, lifters_capacities):
    # Try all possible combinations of lifters
    for r in range(1, len(lifters_capacities) + 1):
        for combo in itertools.combinations(range(len(lifters_capacities)), r):
            if sum(lifters_capacities[i] for i in combo) >= box_weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [44, 70, 139, 184, 194, 67, 188, 33, 149, 194, 132, 117]
    lifters_capacities = [110, 109, 30, 34]
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters_capacities)))
        
        # Sort boxes in descending order for better optimization
        boxes_to_process = sorted([(w, i) for i, w in enumerate(remaining_boxes)], reverse=True)
        
        for box_weight, original_index in boxes_to_process:
            if not available_lifters:
                break
                
            # Try to find lifters for current box
            possible_combinations = []
            for r in range(1, len(available_lifters) + 1):
                for combo in itertools.combinations(available_lifters, r):
                    if sum(lifters_capacities[i] for i in combo) >= box_weight:
                        possible_combinations.append(list(combo))
            
            if possible_combinations:
                # Use the combination with minimum number of lifters
                best_combo = min(possible_combinations, key=len)
                step.append((box_weight, sorted(best_combo)))
                available_lifters -= set(best_combo)
                remaining_boxes.remove(box_weight)
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(f"<<{result.strip()}>>")

solve_box_lifting()