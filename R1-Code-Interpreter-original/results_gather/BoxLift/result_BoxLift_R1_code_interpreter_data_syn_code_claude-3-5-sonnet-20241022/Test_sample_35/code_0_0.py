from itertools import combinations

def can_lift(lifters, weight):
    # Check all possible combinations of lifters to see if they can lift the weight
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [16, 56, 15, 82, 85, 12, 12, 32]
    lifters = [63, 48, 72, 68]
    
    # Sort boxes in descending order to handle heavy boxes first
    box_indices = list(range(len(boxes)))
    box_indices.sort(key=lambda x: boxes[x], reverse=True)
    
    solution = []
    remaining_boxes = set(box_indices)
    step = 1
    
    while remaining_boxes and step <= 3:
        step_solution = []
        used_lifters = set()
        
        # Try to assign boxes to lifters
        for box_idx in list(remaining_boxes):
            if used_lifters == set(range(len(lifters))):
                break
                
            # Get available lifters
            available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
            available_lifter_weights = [lifters[i] for i in available_lifters]
            
            # Try to lift the current box
            lifter_combo = can_lift(available_lifter_weights, boxes[box_idx])
            if lifter_combo:
                # Convert local lifter indices to global indices
                global_lifter_indices = [available_lifters[i] for i in lifter_combo]
                step_solution.append((boxes[box_idx], global_lifter_indices))
                remaining_boxes.remove(box_idx)
                used_lifters.update(global_lifter_indices)
        
        if step_solution:
            solution.append(f"Step {step}: {step_solution}")
        step += 1
    
    if remaining_boxes:
        return "No solution found within 3 steps"
    return "\n".join(solution)

print(solve_box_lifting())