from itertools import combinations
import copy

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(weight, available_lifters, max_lifters=6):
    lifter_capacities = [available_lifters[i] for i in range(len(available_lifters))]
    valid_combinations = []
    
    # Try combinations of different sizes
    for size in range(1, min(len(lifter_capacities) + 1, max_lifters + 1)):
        for combo_indices in combinations(range(len(lifter_capacities)), size):
            combo_capacities = [lifter_capacities[i] for i in combo_indices]
            if sum(combo_capacities) >= weight:
                valid_combinations.append(list(combo_indices))
    
    # Sort by number of lifters (prefer fewer lifters)
    valid_combinations.sort(key=len)
    return valid_combinations

def solve_box_lifting():
    boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
    lifters = [76, 78, 96, 122, 74, 80]
    
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    steps = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box_idx, box_weight in remaining_boxes[:]:
            # Find valid lifter combinations for this box
            valid_combinations = find_lifter_combinations(box_weight, 
                                                        [lifters[i] for i in available_lifters])
            
            # Try each combination
            assigned = False
            for combo in valid_combinations:
                # Check if all lifters in combo are still available
                if all(l in available_lifters for l in combo):
                    # Assign the box to these lifters
                    step.append((box_weight, [available_lifters[i] for i in combo]))
                    # Remove used lifters
                    for l in combo:
                        available_lifters.remove(available_lifters[l])
                    remaining_boxes.remove((box_idx, box_weight))
                    assigned = True
                    break
                    
            if not available_lifters:
                break
                
        steps.append(step)
        
        if len(steps) > 7:
            return None
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())