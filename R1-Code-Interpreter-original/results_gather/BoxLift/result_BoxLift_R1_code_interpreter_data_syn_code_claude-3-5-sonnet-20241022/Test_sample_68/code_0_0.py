from itertools import combinations
import copy

def can_lift(lifters_indices, lifters_capacity, weight):
    return sum(lifters_capacity[i] for i in lifters_indices) >= weight

def find_solution():
    boxes = [190, 107, 357, 385, 343, 392, 97, 277, 231, 77, 147, 138, 371, 303, 173, 397, 123, 72, 241, 291]
    lifters = [100, 70, 149, 198, 118, 100, 153]
    
    # Sort boxes in descending order with their original indices
    boxes_with_idx = list(enumerate(boxes))
    boxes_with_idx.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_idx.copy()
    
    while remaining_boxes:
        step_assignments = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            assigned = False
            for size in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifter_combo in combinations(available_lifters, size):
                    if can_lift(lifter_combo, lifters, box_weight):
                        step_assignments.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.remove((box_idx, box_weight))
                        assigned = True
                        break
        
        if step_assignments:
            solution.append(step_assignments)
        
        if len(solution) > 7:
            return None
    
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

solution = find_solution()
if solution:
    print(solution)
else:
    print("No solution found within 7 steps")