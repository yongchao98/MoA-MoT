from itertools import combinations

def can_lift_box(lifters_capacity, box_weight):
    return sum(lifters_capacity) >= box_weight

def find_best_assignments(boxes, lifters, used_lifters=None):
    if used_lifters is None:
        used_lifters = set()
    
    assignments = []
    remaining_boxes = boxes.copy()
    available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
    
    # Try single lifters first
    for lifter_idx in available_lifters:
        for box in remaining_boxes:
            if lifters[lifter_idx] >= box:
                assignments.append((box, [lifter_idx]))
                remaining_boxes.remove(box)
                break
    
    # Try combinations for heavier boxes
    for box in remaining_boxes.copy():
        for r in range(2, len(available_lifters) + 1):
            for lifter_combo in combinations(available_lifters, r):
                total_capacity = sum(lifters[i] for i in lifter_combo)
                if total_capacity >= box:
                    assignments.append((box, list(lifter_combo)))
                    remaining_boxes.remove(box)
                    break
            if box not in remaining_boxes:
                break
                
    return assignments

def solve_box_lifting():
    boxes = [280, 32, 168, 120, 112, 285, 46, 33, 148, 74, 295, 279, 207, 138, 153, 63]
    lifters = [53, 92, 104, 86, 106, 68]
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 6:
        used_lifters = set()
        step_assignments = []
        
        # Sort boxes in descending order
        remaining_boxes.sort(reverse=True)
        
        while remaining_boxes:
            possible_assignments = find_best_assignments([b for b in remaining_boxes], 
                                                       lifters, used_lifters)
            if not possible_assignments:
                break
            
            for box_weight, lifter_indices in possible_assignments:
                if not any(idx in used_lifters for idx in lifter_indices):
                    step_assignments.append((box_weight, lifter_indices))
                    used_lifters.update(lifter_indices)
                    remaining_boxes.remove(box_weight)
                    
                if len(used_lifters) == len(lifters):
                    break
            
            if len(used_lifters) == len(lifters):
                break
        
        if step_assignments:
            solution.append(step_assignments)
    
    # Format and print solution
    result = ""
    for i, step_assignments in enumerate(solution, 1):
        result += f"Step {i}: {step_assignments}\n"
    print(result.strip())

solve_box_lifting()