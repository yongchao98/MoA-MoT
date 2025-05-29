from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Try all possible combinations of lifters
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [17, 69, 16, 83, 95, 91, 68, 25]
    lifter_caps = [42, 68, 50, 48]
    
    # Sort boxes in descending order for better efficiency
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index[:]
    step = 1
    
    while remaining_boxes and step <= 4:
        step_solution = []
        available_lifters = set(range(len(lifter_caps)))
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            box_idx, weight = remaining_boxes[i]
            
            # Find lifters that can handle this box
            possible_lifters = []
            for r in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, r):
                    if sum(lifter_caps[j] for j in combo) >= weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step_solution.append((weight, list(possible_lifters)))
                available_lifters -= set(possible_lifters)
                remaining_boxes.pop(i)
            else:
                i += 1
        
        if step_solution:
            solution.append(f"Step {step}: {step_solution}")
        step += 1
    
    if remaining_boxes:
        return "No solution found"
    return "\n".join(solution)

print(solve_box_lifting())