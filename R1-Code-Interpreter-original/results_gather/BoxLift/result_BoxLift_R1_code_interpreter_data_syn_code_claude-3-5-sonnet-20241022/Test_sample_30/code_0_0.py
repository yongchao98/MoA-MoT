from itertools import combinations

def can_lift_box(box_weight, lifters_capacities):
    # Returns all possible combinations of lifters that can lift this box
    all_possible = []
    for i in range(1, len(lifters_capacities) + 1):
        for combo in combinations(range(len(lifters_capacities)), i):
            if sum(lifters_capacities[j] for j in combo) >= box_weight:
                all_possible.append(list(combo))
    return all_possible

def find_solution():
    boxes = [14, 51, 14, 69, 77, 20, 52, 96]
    lifters = [69, 66, 74]
    box_indices = list(range(len(boxes)))
    solution = []
    
    while boxes:
        step = []
        used_lifters = set()
        remaining_boxes = boxes.copy()
        remaining_indices = box_indices.copy()
        
        # Try to assign boxes to lifters in current step
        while remaining_boxes and len(used_lifters) < len(lifters):
            # Find heaviest box that can be lifted with remaining lifters
            best_box = None
            best_box_idx = None
            best_lifters = None
            
            for box, box_idx in zip(remaining_boxes, remaining_indices):
                available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
                available_capacities = [lifters[i] for i in available_lifters]
                possible_lifters = can_lift_box(box, available_capacities)
                
                if possible_lifters:
                    # Convert local lifter indices to global indices
                    global_lifters = [[available_lifters[i] for i in combo] for combo in possible_lifters]
                    if best_box is None or box > best_box:
                        best_box = box
                        best_box_idx = box_idx
                        best_lifters = global_lifters[0]  # Take first valid combination
            
            if best_box is not None:
                step.append((best_box, best_lifters))
                used_lifters.update(best_lifters)
                idx = remaining_boxes.index(best_box)
                remaining_boxes.pop(idx)
                remaining_indices.pop(idx)
                boxes.remove(best_box)
                box_indices.remove(best_box_idx)
            else:
                break
        
        if step:
            solution.append(step)
        else:
            return None
    
    return solution

solution = find_solution()
if solution:
    for i, step in enumerate(solution, 1):
        print(f"Step {i}: {step}")