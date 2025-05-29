def find_lifting_solution():
    boxes = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
    lifters = [184, 149, 64, 190, 59, 51]
    
    # Sort boxes in descending order and keep track of original indices
    boxes = [(i, w) for i, w in enumerate(boxes)]
    boxes.sort(key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    solution = []
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifters)))
        
        # Process remaining boxes
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            box_idx, weight = remaining_boxes[i]
            
            # Try combinations of lifters for current box
            for r in range(1, len(available_lifters) + 1):
                lifter_combinations = itertools.combinations(available_lifters, r)
                assigned = False
                
                for lifters_combo in lifter_combinations:
                    if sum(lifters[j] for j in lifters_combo) >= weight:
                        step.append((weight, list(lifters_combo)))
                        for j in lifters_combo:
                            available_lifters.remove(j)
                        remaining_boxes.pop(i)
                        assigned = True
                        break
                
                if assigned:
                    break
            
            if not assigned:
                i += 1
        
        if step:
            solution.append(step)
        
        if len(solution) > 8:
            return "No solution found within 8 steps"
    
    # Format solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

print(find_lifting_solution())