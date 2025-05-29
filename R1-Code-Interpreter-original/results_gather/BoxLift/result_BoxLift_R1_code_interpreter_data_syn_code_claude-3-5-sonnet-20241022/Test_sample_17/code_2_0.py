from itertools import combinations

def solve_box_lifting():
    boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
    lifters = [76, 78, 96, 122, 74, 80]
    
    def can_lift(lifter_indices, weight):
        return sum(lifters[i] for i in lifter_indices) >= weight
    
    # Sort boxes by weight (descending)
    remaining_boxes = [(i, w) for i, w in enumerate(boxes)]
    remaining_boxes.sort(key=lambda x: x[1], reverse=True)
    
    steps = []
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        i = 0
        while i < len(remaining_boxes):
            _, weight = remaining_boxes[i]
            assigned = False
            
            # Try different combinations of available lifters
            for num_lifters in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifter_combo in combinations(available_lifters, num_lifters):
                    if can_lift(lifter_combo, weight):
                        step.append((weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
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
    result = []
    for i, step in enumerate(steps, 1):
        result.append(f"Step {i}: {step}")
    return "\n".join(result)

print(solve_box_lifting())