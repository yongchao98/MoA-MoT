import itertools

# Input data
boxes = [356, 292, 240, 155, 142, 223, 178, 231, 254, 315, 236, 328, 203, 397, 44, 200, 80, 131, 175, 216]
lifters = [155, 116, 53, 98, 96, 156, 167]

def can_lift(lifters_subset, weight):
    return sum(lifters[i] for i in lifters_subset) >= weight

def find_solution():
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    remaining_boxes = boxes_with_index.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters in this step
        while available_lifters and remaining_boxes:
            found = False
            current_box = None
            best_lifters = None
            
            # Try each remaining box
            for box_idx, (box_orig_idx, box_weight) in enumerate(remaining_boxes):
                # Try different combinations of available lifters
                for r in range(1, len(available_lifters) + 1):
                    for lifter_combination in itertools.combinations(available_lifters, r):
                        if can_lift(lifter_combination, box_weight):
                            current_box = remaining_boxes.pop(box_idx)
                            best_lifters = list(lifter_combination)
                            available_lifters -= set(lifter_combination)
                            step.append((box_weight, best_lifters))
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            
            if not found:
                break
                
        steps.append(step)
        if len(steps) > 7:
            return None
    
    return steps

# Get and format solution
solution = find_solution()
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(output.strip())
else:
    print("No solution found within 7 steps")