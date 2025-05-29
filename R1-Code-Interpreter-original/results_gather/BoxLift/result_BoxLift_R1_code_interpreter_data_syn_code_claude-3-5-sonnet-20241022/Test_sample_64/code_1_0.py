import itertools

def solve_box_lifting():
    boxes = [45, 139, 291, 189, 171, 265, 139, 312, 289, 394, 183, 229, 321, 250, 226, 335, 59, 127, 298, 339, 256, 161, 288, 179]
    lifters = [52, 99, 86, 56, 92, 134]
    
    # Sort boxes from heaviest to lightest
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    remaining_boxes = boxes_with_index.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters, prioritizing efficient combinations
        remaining_boxes_copy = remaining_boxes.copy()
        
        while remaining_boxes_copy and available_lifters:
            best_assignment = None
            best_efficiency = 0
            
            for box_idx, box_weight in remaining_boxes_copy:
                for r in range(1, len(available_lifters) + 1):
                    for combo in itertools.combinations(available_lifters, r):
                        total_capacity = sum(lifters[i] for i in combo)
                        if total_capacity >= box_weight:
                            efficiency = box_weight / total_capacity
                            if efficiency > best_efficiency:
                                best_efficiency = efficiency
                                best_assignment = (box_idx, box_weight, list(combo))
                            break
                    if best_efficiency > 0.8:  # Good enough efficiency threshold
                        break
            
            if best_assignment:
                box_idx, box_weight, lifter_combo = best_assignment
                step.append((box_weight, lifter_combo))
                available_lifters -= set(lifter_combo)
                remaining_boxes_copy = [(i, w) for i, w in remaining_boxes_copy if i != box_idx]
                remaining_boxes = [(i, w) for i, w in remaining_boxes if i != box_idx]
            else:
                break
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format and print the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()