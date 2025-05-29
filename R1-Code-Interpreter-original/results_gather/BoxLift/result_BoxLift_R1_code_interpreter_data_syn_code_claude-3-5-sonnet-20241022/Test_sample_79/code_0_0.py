from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Returns list of lifters that can lift the weight together
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [368, 220, 131, 287, 153, 175, 153, 217, 262, 115, 365, 394, 
             291, 222, 142, 126, 189, 82, 159, 367, 356, 168, 188, 355]
    lifter_capacities = [200, 91, 128, 198, 196, 119, 188]
    
    # Sort boxes by weight (descending)
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find lifters that can handle this box
            possible_lifters = []
            for r in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, r):
                    if sum(lifter_capacities[i] for i in combo) >= box_weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step.append((box_weight, possible_lifters))
                available_lifters -= set(possible_lifters)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        
        if len(solution) > 6:
            return None
    
    # Format solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()