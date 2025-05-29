from itertools import combinations
import copy

def can_lift(lifters, weight):
    if not lifters:
        return False
    return sum(lifters) >= weight

def find_minimal_lifters(lifter_capacities, weight, used_lifters):
    available_lifters = [i for i in range(len(lifter_capacities)) if i not in used_lifters]
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(lifter_capacities[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [75, 210, 297, 242, 276, 134, 383, 358, 146, 394, 59, 350, 94, 152, 250, 136, 116, 342, 283, 87]
    lifters = [101, 62, 125, 126, 138, 54]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign boxes in current step
        remaining_boxes_copy = remaining_boxes.copy()
        for box_idx, box_weight in remaining_boxes_copy:
            lifters_for_box = find_minimal_lifters(lifters, box_weight, used_lifters)
            if lifters_for_box:
                step.append((box_weight, lifters_for_box))
                used_lifters.update(lifters_for_box)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            solution.append(step)
        
        if len(solution) > 10:
            return "No solution found within 10 steps"
    
    # Format the solution
    formatted_solution = ""
    for i, step in enumerate(solution, 1):
        formatted_solution += f"Step {i}: {step}\n"
    
    return formatted_solution.strip()

print(solve_box_lifting())