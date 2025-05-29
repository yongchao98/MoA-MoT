def solve_box_lifting():
    boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
    lifters = [76, 78, 96, 122, 74, 80]
    
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    def get_lifter_combinations(weight, available_lifters):
        results = []
        for i in range(1, len(available_lifters) + 1):
            for combo in combinations(available_lifters, i):
                if sum(lifters[j] for j in combo) >= weight:
                    results.append(list(combo))
        return sorted(results, key=len)  # Prefer fewer lifters
    
    steps = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        i = 0
        while i < len(remaining_boxes):
            box_idx, box_weight = remaining_boxes[i]
            valid_combos = get_lifter_combinations(box_weight, list(available_lifters))
            
            assigned = False
            for combo in valid_combos:
                if all(l in available_lifters for l in combo):
                    step.append((box_weight, combo))
                    available_lifters -= set(combo)
                    remaining_boxes.pop(i)
                    assigned = True
                    break
            
            if not assigned:
                i += 1
                
            if not available_lifters:
                break
        
        if step:
            steps.append(step)
        
        if len(steps) > 7:
            return "No solution found within 7 steps"
        
        if not step and remaining_boxes:
            return "Failed to assign all boxes"
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

print(solve_box_lifting())