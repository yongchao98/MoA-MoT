import itertools

def find_lifter_combinations(lifters, target_weight):
    # Find all possible combinations of lifters that can lift a given weight
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for comb in itertools.combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in comb) >= target_weight:
                valid_combinations.append(list(comb))
    return sorted(valid_combinations, key=len)  # Sort by number of lifters needed

def solve_box_lifting():
    boxes = [260, 319, 78, 313, 164, 297, 367, 392, 175, 245, 200, 45, 386, 124, 169, 97, 359, 273, 175, 238]
    lifters = [57, 82, 155, 87, 174, 192, 149]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    solution = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combs = [c for c in find_lifter_combinations(lifters, box) 
                         if all(i in available_lifters for i in c)]
            
            if valid_combs:
                best_comb = valid_combs[0]  # Take the combination with fewest lifters
                step.append((box, best_comb))
                remaining_boxes.remove(box)
                available_lifters -= set(best_comb)
        
        if step:
            solution.append(step)
    
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())