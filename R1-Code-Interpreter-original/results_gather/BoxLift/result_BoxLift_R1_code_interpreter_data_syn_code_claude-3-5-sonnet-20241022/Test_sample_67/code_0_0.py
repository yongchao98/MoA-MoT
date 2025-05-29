from itertools import combinations
import copy

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(lifters, weight, max_lifters=6):
    valid_combinations = []
    for i in range(1, min(len(lifters) + 1, max_lifters + 1)):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                valid_combinations.append(list(combo))
    return sorted(valid_combinations, key=len)

def solve_box_lifting():
    boxes = [93, 50, 31, 217, 183, 34, 43, 268, 281, 93, 145, 74, 278, 272, 86, 81]
    lifters = [97, 53, 143, 47, 94, 132]
    
    # Sort boxes from heaviest to lightest
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    steps = []
    remaining_boxes = boxes_with_index[:]
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combos = [c for c in find_lifter_combinations(lifters, box_weight)
                          if all(l in available_lifters for l in c)]
            
            if valid_combos:
                best_combo = min(valid_combos, key=len)
                step.append((box_weight, best_combo))
                available_lifters -= set(best_combo)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()