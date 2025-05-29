from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Returns list of lifter combinations that can lift the weight
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for comb in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in comb) >= weight:
                valid_combinations.append(list(comb))
    return valid_combinations

def solve_box_lifting():
    boxes = [94, 177, 88, 131, 87, 113, 148, 86, 141, 83, 86, 73]
    lifters = [55, 92, 35, 48, 96]
    
    # Sort boxes by weight (descending)
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combs = [c for c in can_lift(lifters, box_weight) 
                         if all(l in available_lifters for l in c)]
            
            if valid_combs:
                # Use the combination with minimum number of lifters
                best_comb = min(valid_combs, key=len)
                step.append((box_weight, best_comb))
                available_lifters -= set(best_comb)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        
        if len(solution) > 6:
            return None
    
    # Format solution
    result = []
    for i, step in enumerate(solution, 1):
        result.append(f"Step {i}: {step}")
    
    return "\n".join(result)

print(solve_box_lifting())