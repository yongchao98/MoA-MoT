from itertools import combinations

def get_valid_lifter_combinations(weight, lifters):
    valid_combos = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                valid_combos.append(list(combo))
    return valid_combos

def solve_box_lifting():
    boxes = [128, 144, 42, 49, 132, 118, 32, 48, 21, 94, 85, 51]
    lifter_capacities = [85, 90, 40, 102]
    
    # Sort boxes from heaviest to lightest
    sorted_boxes = sorted(enumerate(boxes), key=lambda x: -x[1])
    solution = []
    remaining_boxes = sorted_boxes.copy()
    
    while remaining_boxes:
        current_step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # First handle heavy boxes that require multiple lifters
        for box_idx, weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            valid_combinations = get_valid_lifter_combinations(
                weight,
                [lifter_capacities[i] for i in available_lifters]
            )
            
            if valid_combinations:
                # Choose the combination that uses minimum lifters
                best_combo = min(valid_combinations, key=len)
                if all(l in available_lifters for l in best_combo):
                    current_step.append((weight, best_combo))
                    for lifter in best_combo:
                        available_lifters.remove(lifter)
                    remaining_boxes.remove((box_idx, weight))
        
        if current_step:
            solution.append(current_step)
        
        # Break if we exceed 5 steps
        if len(solution) > 5:
            print("<<<No solution found within 5 steps>>>")
            return
        
        # Break if we can't make progress
        if not current_step and remaining_boxes:
            print("<<<No solution found within 5 steps>>>")
            return
    
    # Format and print solution
    result = []
    for i, step in enumerate(solution, 1):
        result.append(f"Step {i}: {step}")
    print("<<<" + "\n".join(result) + ">>>")

solve_box_lifting()