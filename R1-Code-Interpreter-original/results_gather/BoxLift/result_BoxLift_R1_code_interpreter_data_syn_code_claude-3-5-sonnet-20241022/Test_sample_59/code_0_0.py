import copy

def find_lifter_combinations(lifters, target_weight, max_lifters=7):
    n = len(lifters)
    best_combination = None
    min_lifters_used = float('inf')
    
    # Try combinations of lifters
    for i in range(1, min(max_lifters + 1, n + 1)):
        def backtrack(start, current_sum, current_combination):
            nonlocal best_combination, min_lifters_used
            
            if current_sum >= target_weight and len(current_combination) < min_lifters_used:
                best_combination = current_combination[:]
                min_lifters_used = len(current_combination)
                return
            
            if len(current_combination) >= i:
                return
                
            for j in range(start, n):
                current_combination.append(j)
                backtrack(j + 1, current_sum + lifters[j], current_combination)
                current_combination.pop()
                
        backtrack(0, 0, [])
        if best_combination is not None:
            return best_combination
            
    return None

def solve_box_lifting(boxes, lifters):
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    lifter_capacities = lifters
    remaining_boxes = boxes[:]
    steps = []
    
    while remaining_boxes and len(steps) < 8:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        current_boxes = remaining_boxes[:]
        for box_weight, box_idx in current_boxes:
            if not available_lifters:
                break
                
            # Find lifters that can lift this box
            lifters_needed = find_lifter_combinations(
                [lifter_capacities[i] for i in available_lifters],
                box_weight,
                len(available_lifters)
            )
            
            if lifters_needed is not None:
                # Convert relative indices to absolute indices
                absolute_lifters = [list(available_lifters)[i] for i in lifters_needed]
                step.append((box_weight, absolute_lifters))
                available_lifters -= set(absolute_lifters)
                remaining_boxes.remove((box_weight, box_idx))
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result if len(remaining_boxes) == 0 else "No solution found"

# Input data
boxes = [299, 338, 270, 207, 76, 195, 101, 89, 277, 236, 240, 342, 47, 207, 220, 349, 185, 124, 309, 112, 266, 278, 162, 119]
lifters = [198, 55, 76, 137, 163, 117, 146]

print(solve_box_lifting(boxes, lifters))