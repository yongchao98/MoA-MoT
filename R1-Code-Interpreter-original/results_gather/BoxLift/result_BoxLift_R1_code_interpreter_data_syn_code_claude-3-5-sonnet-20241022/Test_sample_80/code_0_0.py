from itertools import combinations

def can_lift(lifters, box_weight):
    # Check if given lifters can lift the box
    return sum(lifters) >= box_weight

def find_solution():
    boxes = [75, 73, 64, 59, 58, 47, 14, 13]  # Sorted in descending order
    lifter_capacities = [74, 57, 51]
    n_lifters = len(lifter_capacities)
    
    # Store the solution
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(n_lifters))
        
        # Try to assign boxes to lifters in current step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try single lifters first
            assigned = False
            for lifter in list(available_lifters):
                if lifter_capacities[lifter] >= box:
                    step.append((box, [lifter]))
                    available_lifters.remove(lifter)
                    remaining_boxes.remove(box)
                    assigned = True
                    break
            
            if assigned:
                continue
                
            # Try combinations of available lifters
            for r in range(2, len(available_lifters) + 1):
                if assigned:
                    break
                for lifter_combo in combinations(available_lifters, r):
                    total_capacity = sum(lifter_capacities[i] for i in lifter_combo)
                    if total_capacity >= box:
                        step.append((box, list(lifter_combo)))
                        for lifter in lifter_combo:
                            available_lifters.remove(lifter)
                        remaining_boxes.remove(box)
                        assigned = True
                        break
        
        if step:
            solution.append(step)
    
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(find_solution())