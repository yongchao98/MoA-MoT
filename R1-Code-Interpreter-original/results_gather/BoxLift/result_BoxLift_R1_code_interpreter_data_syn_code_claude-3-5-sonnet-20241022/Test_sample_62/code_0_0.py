from itertools import combinations
import copy

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_solution():
    boxes = [251, 221, 209, 207, 194, 184, 176, 159, 152, 152, 141, 128, 91, 57, 54, 47]
    lifter_capacities = [149, 131, 113, 109, 124]
    n_lifters = len(lifter_capacities)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(n_lifters))
        
        # Try to assign boxes to lifters in this step
        while available_lifters and remaining_boxes:
            assigned = False
            
            # Try single lifter first
            for box_idx, box in enumerate(remaining_boxes):
                for lifter in list(available_lifters):
                    if lifter_capacities[lifter] >= box:
                        step.append((box, [lifter]))
                        available_lifters.remove(lifter)
                        remaining_boxes.pop(box_idx)
                        assigned = True
                        break
                if assigned:
                    break
            
            if not assigned and available_lifters:
                # Try combinations of lifters for heavy boxes
                for n_combined in range(2, len(available_lifters) + 1):
                    if assigned:
                        break
                    for lifter_combo in combinations(available_lifters, n_combined):
                        combined_capacity = sum(lifter_capacities[i] for i in lifter_combo)
                        for box_idx, box in enumerate(remaining_boxes):
                            if combined_capacity >= box:
                                step.append((box, list(lifter_combo)))
                                for lifter in lifter_combo:
                                    available_lifters.remove(lifter)
                                remaining_boxes.pop(box_idx)
                                assigned = True
                                break
                        if assigned:
                            break
            
            if not assigned:
                break
                
        if step:
            solution.append(step)
        
        if len(solution) > 6:
            return None
            
    return solution

solution = find_solution()
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(f"<<{output}>>")
else:
    print("No solution found within 6 steps")