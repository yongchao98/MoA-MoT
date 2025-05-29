from itertools import combinations
import copy

def can_lift_box(lifters_indices, lifters_capacity, box_weight):
    return sum(lifters_capacity[i] for i in lifters_indices) >= box_weight

def find_solution():
    boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
    lifters = [93, 71, 82, 91]
    
    # Sort boxes from heaviest to lightest
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes and len(solution) < 5:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                lifter_combinations = combinations(available_lifters, r)
                
                for lifter_combo in lifter_combinations:
                    if can_lift_box(lifter_combo, lifters, box_weight):
                        step.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.remove((box_idx, box_weight))
                        break
                if (box_weight, list(lifter_combo)) in step:
                    break
        
        if step:
            solution.append(step)
        else:
            break
    
    # Format and print solution
    if not remaining_boxes and len(solution) <= 5:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
        return True
    return False

find_solution()