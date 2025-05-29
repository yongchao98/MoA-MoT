import itertools

def can_lift(weight, lifters_capacity):
    return sum(lifters_capacity) >= weight

def find_lifter_combination(weight, lifters, used_lifters):
    available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
    for r in range(1, len(available_lifters) + 1):
        for comb in itertools.combinations(available_lifters, r):
            if sum(lifters[i] for i in comb) >= weight:
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [75, 298, 70, 200, 278, 35, 278, 56, 270, 212, 161, 98, 244, 292, 177, 145]
    lifters = [105, 102, 125, 93, 123]
    
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes):
            box_idx, weight = remaining_boxes[i]
            lifter_combination = find_lifter_combination(weight, lifters, used_lifters)
            
            if lifter_combination is not None:
                step.append((weight, lifter_combination))
                used_lifters.update(lifter_combination)
                remaining_boxes.pop(i)
            else:
                i += 1
                
            if len(used_lifters) == len(lifters):
                break
        
        solution.append(step)
    
    # Print solution
    for i, step in enumerate(solution, 1):
        print(f"Step {i}: {step}")
    print(f"\nTotal steps: {len(solution)}")

solve_box_lifting()