from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [18, 27, 65, 31, 60, 73, 30, 67]
    lifter_capacities = [66, 72, 43]
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 4:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort remaining boxes in descending order
        boxes_with_indices = list(enumerate(remaining_boxes))
        boxes_with_indices.sort(key=lambda x: x[1], reverse=True)
        
        new_remaining = remaining_boxes.copy()
        
        for orig_idx, weight in boxes_with_indices:
            if not available_lifters:
                break
                
            # Try to lift current box with available lifters
            current_lifters = [i for i in range(len(lifter_capacities)) if i in available_lifters]
            lifter_combo = can_lift([lifter_capacities[i] for i in current_lifters], weight)
            
            if lifter_combo is not None:
                # Convert local lifter indices to global indices
                global_lifter_indices = [current_lifters[i] for i in lifter_combo]
                step.append((weight, global_lifter_indices))
                new_remaining[orig_idx] = 0  # Mark as lifted
                # Remove used lifters
                for lifter in global_lifter_indices:
                    available_lifters.remove(lifter)
        
        if step:
            solution.append(step)
            remaining_boxes = [x for x in new_remaining if x != 0]
        else:
            break
    
    # Format the solution
    if not remaining_boxes and len(solution) <= 4:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        return result.strip()
    return "No solution found within 4 steps"

print(solve_box_lifting())