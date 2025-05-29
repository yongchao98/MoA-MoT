import itertools

def find_lifting_solution():
    boxes = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
    lifters = [184, 149, 64, 190, 59, 51]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    solution = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # First, handle heavy boxes that need multiple lifters
        for box_idx, (orig_idx, weight) in enumerate(remaining_boxes):
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                if weight <= sum(lifters[i] for i in list(available_lifters)[:r]):
                    lifter_combinations = itertools.combinations(available_lifters, r)
                    for lifter_combo in lifter_combinations:
                        if sum(lifters[i] for i in lifter_combo) >= weight:
                            step.append((weight, list(lifter_combo)))
                            available_lifters -= set(lifter_combo)
                            remaining_boxes[box_idx] = None
                            break
                    break
        
        # Remove processed boxes
        remaining_boxes = [box for box in remaining_boxes if box is not None]
        
        # Handle lighter boxes with single lifters
        if available_lifters:
            for box_idx, (orig_idx, weight) in enumerate(remaining_boxes):
                if not available_lifters:
                    break
                    
                for lifter in list(available_lifters):
                    if lifters[lifter] >= weight:
                        step.append((weight, [lifter]))
                        available_lifters.remove(lifter)
                        remaining_boxes[box_idx] = None
                        break
                        
            remaining_boxes = [box for box in remaining_boxes if box is not None]
        
        if step:
            solution.append(step)
            
        if len(solution) > 8:
            return None
    
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(find_lifting_solution())