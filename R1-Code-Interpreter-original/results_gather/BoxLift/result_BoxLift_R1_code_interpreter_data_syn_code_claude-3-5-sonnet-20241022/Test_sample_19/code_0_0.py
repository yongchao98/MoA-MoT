from itertools import combinations
import copy

def can_lift_box(lifters_capacities, box_weight):
    # Try combinations of lifters to lift the box
    for i in range(1, len(lifters_capacities) + 1):
        for combo in combinations(range(len(lifters_capacities)), i):
            if sum(lifters_capacities[j] for j in combo) >= box_weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [55, 77, 274, 297, 70, 312, 108, 84, 340, 273, 157, 386, 95, 223, 311, 306, 287, 91, 220, 296, 320, 237, 97, 73]
    lifters = [167, 138, 133, 87, 142, 167, 140]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find lifters that can lift this box
            lifters_subset = [lifters[i] for i in available_lifters]
            lifter_indices = can_lift_box(lifters_subset, box)
            
            if lifter_indices is not None:
                # Convert relative indices to absolute indices
                absolute_indices = [list(available_lifters)[i] for i in lifter_indices]
                step.append((box, absolute_indices))
                # Remove used lifters from available set
                for idx in absolute_indices:
                    available_lifters.remove(idx)
                remaining_boxes.remove(box)
        
        if step:
            solution.append(step)
        
        if len(solution) > 7:
            return None
    
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

# Solve and print the solution
print(solve_box_lifting())