from itertools import combinations

def can_lift_box(box_weight, lifters_capacities):
    # Try all possible combinations of lifters
    for i in range(1, len(lifters_capacities) + 1):
        for lifter_combo in combinations(range(len(lifters_capacities)), i):
            if sum(lifters_capacities[j] for j in lifter_combo) >= box_weight:
                return list(lifter_combo)
    return None

def find_solution():
    boxes = [41, 31, 60, 20, 11, 26, 52, 98]
    lifter_capacities = [49, 55, 63]
    remaining_boxes = boxes.copy()
    solution = []
    
    while remaining_boxes and len(solution) < 4:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort boxes in descending order for greedy approach
        boxes_to_try = sorted([(w, i) for i, w in enumerate(remaining_boxes)], reverse=True)
        
        for box_weight, original_index in boxes_to_try:
            if not available_lifters:
                break
                
            # Try to lift current box with available lifters
            available_lifter_capacities = [lifter_capacities[i] for i in available_lifters]
            lifters_needed = can_lift_box(box_weight, available_lifter_capacities)
            
            if lifters_needed is not None:
                # Convert local lifter indices to global indices
                global_lifters = [list(available_lifters)[i] for i in lifters_needed]
                step.append((box_weight, global_lifters))
                # Remove used lifters from available set
                for lifter in global_lifters:
                    available_lifters.remove(lifter)
                # Remove box from remaining boxes
                remaining_boxes.remove(box_weight)
        
        if step:
            solution.append(step)
        else:
            break
    
    if not remaining_boxes and len(solution) <= 4:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
        return True
    return False

find_solution()