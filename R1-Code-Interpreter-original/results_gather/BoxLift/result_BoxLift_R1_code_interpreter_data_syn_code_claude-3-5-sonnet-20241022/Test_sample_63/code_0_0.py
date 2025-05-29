import itertools

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(lifter_capacities, weight, max_lifters=6):
    valid_combinations = []
    for i in range(1, max_lifters + 1):
        for combo in itertools.combinations(range(len(lifter_capacities)), i):
            lifters = [lifter_capacities[j] for j in combo]
            if can_lift(lifters, weight):
                valid_combinations.append(list(combo))
    return valid_combinations

def solve_box_lifting():
    boxes = [49, 186, 267, 243, 352, 74, 160, 115, 138, 301, 250, 145, 294, 232, 144, 293, 287, 358, 267, 266]
    lifter_capacities = [60, 92, 94, 68, 72, 79]
    
    # Sort boxes by weight (descending)
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign as many boxes as possible in this step
        for weight, box_idx in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find all possible lifter combinations for this box
            valid_combos = [c for c in find_lifter_combinations(lifter_capacities, weight)
                          if all(l in available_lifters for l in c)]
            
            if valid_combos:
                # Use the combination with minimum number of lifters
                best_combo = min(valid_combos, key=len)
                step.append((weight, best_combo))
                available_lifters -= set(best_combo)
                remaining_boxes.remove((weight, box_idx))
        
        if step:
            solution.append(step)
        
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()