def solve_box_lifting():
    boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
    capacities = [122, 116, 51, 165, 95, 150, 179]
    
    # Pre-calculate all valid lifter pairs
    valid_pairs = {}
    for weight in sorted(set(boxes), reverse=True):
        valid_pairs[weight] = []
        # Single lifter
        for i in range(len(capacities)):
            if capacities[i] >= weight:
                valid_pairs[weight].append([i])
        # Lifter pairs
        for i in range(len(capacities)):
            for j in range(i+1, len(capacities)):
                if capacities[i] + capacities[j] >= weight:
                    valid_pairs[weight].append([i, j])
    
    # Sort boxes by weight
    box_indices = list(range(len(boxes)))
    box_indices.sort(key=lambda x: boxes[x], reverse=True)
    
    solution = []
    unassigned = set(box_indices)
    
    while unassigned:
        step = []
        available_lifters = set(range(len(capacities)))
        
        # Try to assign boxes in current step
        for box_idx in list(unassigned):
            weight = boxes[box_idx]
            
            # Find valid lifter combination from available lifters
            for lifters in valid_pairs[weight]:
                if all(l in available_lifters for l in lifters):
                    step.append((weight, lifters))
                    unassigned.remove(box_idx)
                    for l in lifters:
                        available_lifters.remove(l)
                    break
            
            if not available_lifters:
                break
        
        if step:
            solution.append(step)
            if len(solution) > 7:
                return None
        else:
            return None
    
    # Format and print solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    print(result.strip())
    return solution

solve_box_lifting()