from itertools import combinations

def can_lift_box(lifters, box_weight):
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for comb in combinations(range(len(lifters)), r):
            total_capacity = sum(lifters[i] for i in comb)
            if total_capacity >= box_weight:
                valid_combinations.append((list(comb), total_capacity))
    return sorted(valid_combinations, key=lambda x: (len(x[0]), -x[1]))

def solve_box_lifting():
    boxes = [196, 90, 142, 22, 143, 101, 43, 166, 79, 129, 77, 187]
    lifters = [71, 58, 104, 105]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 6:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters
        remaining_boxes_copy = remaining_boxes.copy()
        for box_idx, box_weight in remaining_boxes_copy:
            if not available_lifters:
                break
                
            valid_combs = [(c, cap) for c, cap in can_lift_box([lifters[i] for i in available_lifters], box_weight)
                          if all(i in available_lifters for i in c)]
            
            if valid_combs:
                best_comb, capacity = valid_combs[0]
                if capacity >= box_weight:  # Double check capacity
                    step.append((box_weight, best_comb))
                    available_lifters -= set(best_comb)
                    remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        else:
            break
    
    # Verify solution
    for step in solution:
        for box_weight, lifter_indices in step:
            if sum(lifters[i] for i in lifter_indices) < box_weight:
                print("<<<Invalid solution: Lifter capacity exceeded>>>")
                return
    
    # Format and print solution
    if not remaining_boxes and len(solution) <= 6:
        result = []
        for i, step in enumerate(solution, 1):
            step_str = f"Step {i}: {step}"
            result.append(step_str)
        print("<<<" + "\n".join(result) + ">>>")
    else:
        print("<<<No solution found within 6 steps>>>")

solve_box_lifting()