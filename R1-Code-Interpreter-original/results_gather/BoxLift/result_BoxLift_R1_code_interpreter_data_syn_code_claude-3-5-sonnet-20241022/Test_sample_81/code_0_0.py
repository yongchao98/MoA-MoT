from itertools import combinations
import copy

def can_lift_box(lifters_combo, box_weight):
    return sum(lifters_combo) >= box_weight

def find_solution():
    boxes = [174, 117, 74, 120, 98, 65, 112, 61, 46, 87, 176, 110]
    lifters = [119, 104, 67, 98]
    boxes.sort(reverse=True)  # Sort boxes in descending order
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(steps) < 5:
        step = []
        available_lifters = list(range(len(lifters)))
        current_boxes = remaining_boxes.copy()
        
        while current_boxes and available_lifters:
            box = current_boxes[0]
            found = False
            
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                if found:
                    break
                for lifter_combo_idx in combinations(available_lifters, r):
                    lifter_combo = [lifters[i] for i in lifter_combo_idx]
                    if can_lift_box(lifter_combo, box):
                        step.append((box, list(lifter_combo_idx)))
                        for idx in lifter_combo_idx:
                            available_lifters.remove(idx)
                        current_boxes.remove(box)
                        found = True
                        break
            
            if not found:
                break
                
        if step:
            steps.append(step)
            remaining_boxes = current_boxes
        else:
            break
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    if not remaining_boxes and len(steps) <= 5:
        print(result.strip())
        return True
    return False

find_solution()