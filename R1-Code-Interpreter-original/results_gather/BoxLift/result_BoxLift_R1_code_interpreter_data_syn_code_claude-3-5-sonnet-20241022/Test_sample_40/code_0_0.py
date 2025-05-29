from itertools import combinations

def can_lift_box(lifters_combo, box_weight):
    return sum(lifters_combo) >= box_weight

def find_solution():
    boxes = [12, 90, 67, 10, 10, 94, 60, 46]
    lifters = [56, 77, 72]
    boxes_indices = list(range(len(boxes)))
    lifters_indices = list(range(len(lifters)))
    
    # Sort boxes in descending order for better allocation
    sorted_boxes = [(w, i) for i, w in enumerate(boxes)]
    sorted_boxes.sort(reverse=True)
    
    steps = []
    remaining_boxes = set(range(len(boxes)))
    
    while remaining_boxes:
        step = []
        available_lifters = set(lifters_indices)
        
        # Try to assign heaviest boxes first
        for box_weight, box_idx in sorted_boxes:
            if box_idx not in remaining_boxes:
                continue
                
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                lifter_combinations = list(combinations(available_lifters, r))
                for lifter_combo_idx in lifter_combinations:
                    lifter_combo = [lifters[i] for i in lifter_combo_idx]
                    
                    if can_lift_box(lifter_combo, box_weight):
                        step.append((box_weight, list(lifter_combo_idx)))
                        remaining_boxes.remove(box_idx)
                        available_lifters -= set(lifter_combo_idx)
                        break
                if box_idx not in remaining_boxes:
                    break
                    
        steps.append(step)
        if len(steps) > 4:  # Check if we exceed maximum steps
            return None
            
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

solution = find_solution()
print(solution)