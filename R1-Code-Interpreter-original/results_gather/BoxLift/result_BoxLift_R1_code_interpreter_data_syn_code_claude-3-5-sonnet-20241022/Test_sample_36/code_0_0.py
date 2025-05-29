from itertools import combinations

def find_lifter_combinations(box_weight, lifters, used_lifters):
    available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
    valid_combinations = []
    
    # Try all possible combinations of available lifters
    for r in range(1, len(available_lifters) + 1):
        for comb in combinations(available_lifters, r):
            if sum(lifters[i] for i in comb) >= box_weight:
                valid_combinations.append(list(comb))
    
    # Sort by number of lifters used (prefer fewer lifters)
    valid_combinations.sort(key=len)
    return valid_combinations

def solve_box_lifting():
    boxes = [85, 162, 147, 83, 142, 96, 200, 172, 151, 77, 59, 39]
    lifters = [53, 54, 76, 48, 97]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 6:
        step = []
        used_lifters = set()
        
        # Try to assign as many boxes as possible in current step
        i = 0
        while i < len(remaining_boxes):
            box_weight, box_index = remaining_boxes[i]
            combinations = find_lifter_combinations(box_weight, lifters, used_lifters)
            
            if combinations:
                # Use the first valid combination (using minimum lifters)
                lifter_group = combinations[0]
                step.append((box_weight, lifter_group))
                used_lifters.update(lifter_group)
                remaining_boxes.pop(i)
            else:
                i += 1
        
        if step:
            solution.append(step)
        else:
            break
    
    # Format and print solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()