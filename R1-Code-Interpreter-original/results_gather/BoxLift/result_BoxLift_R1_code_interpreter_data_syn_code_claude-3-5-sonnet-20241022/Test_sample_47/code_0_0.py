from itertools import combinations
import copy

def can_lift_box(lifters, box_weight):
    return sum(lifters) >= box_weight

def find_lifter_combinations(lifters, boxes, used_lifters=None):
    if used_lifters is None:
        used_lifters = set()
    
    available_lifters = [cap for i, cap in enumerate(lifters) if i not in used_lifters]
    available_indices = [i for i in range(len(lifters)) if i not in used_lifters]
    assignments = []
    
    # Try to assign single lifters first
    for i, cap in zip(available_indices, available_lifters):
        for box in boxes:
            if cap >= box:
                assignments.append((box, [i]))
                
    # Try combinations for heavier boxes
    for box in boxes:
        for r in range(2, len(available_lifters) + 1):
            for lifter_combo_idx in combinations(range(len(available_lifters)), r):
                lifter_combo = [available_indices[i] for i in lifter_combo_idx]
                total_capacity = sum(lifters[i] for i in lifter_combo)
                if total_capacity >= box:
                    assignments.append((box, lifter_combo))
                    
    return assignments

def solve_box_lifting():
    boxes = [280, 32, 168, 120, 112, 285, 46, 33, 148, 74, 295, 279, 207, 138, 153, 63]
    lifters = [53, 92, 104, 86, 106, 68]
    
    solution = []
    remaining_boxes = boxes.copy()
    step = 1
    
    while remaining_boxes and step <= 6:
        step_assignments = []
        used_lifters = set()
        
        # Sort boxes in descending order to handle heavy boxes first
        remaining_boxes.sort(reverse=True)
        
        while remaining_boxes:
            possible_assignments = find_lifter_combinations(lifters, remaining_boxes, used_lifters)
            if not possible_assignments:
                break
                
            # Choose the assignment that uses the most efficient combination of lifters
            best_assignment = min(possible_assignments, 
                                key=lambda x: (len(x[1]), -x[0]))  # Prefer fewer lifters and heavier boxes
            
            box_weight, lifter_indices = best_assignment
            
            # Update used lifters and remaining boxes
            used_lifters.update(lifter_indices)
            step_assignments.append((box_weight, lifter_indices))
            remaining_boxes.remove(box_weight)
            
            if len(used_lifters) == len(lifters):
                break
        
        solution.append(step_assignments)
        step += 1
    
    # Format and print solution
    result = ""
    for i, step_assignments in enumerate(solution, 1):
        result += f"Step {i}: {step_assignments}\n"
    
    print(result.strip())

solve_box_lifting()