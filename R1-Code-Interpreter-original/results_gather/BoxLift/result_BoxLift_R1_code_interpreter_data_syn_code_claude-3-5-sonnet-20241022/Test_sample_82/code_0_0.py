from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [99, 159, 55, 86, 215, 56, 271, 167, 212, 87, 78, 211, 246, 153, 227, 287]
    lifter_capacities = [143, 50, 70, 79, 136]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    steps = []
    
    while boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters
        remaining_boxes = []
        for weight, box_idx in boxes:
            if not available_lifters:
                remaining_boxes.append((weight, box_idx))
                continue
                
            # Find lifters that can handle this box
            lifters_subset = [lifter_capacities[i] for i in available_lifters]
            lifter_indices = can_lift([lifter_capacities[i] for i in available_lifters], weight)
            
            if lifter_indices is not None:
                # Convert relative indices to absolute indices
                abs_lifter_indices = [list(available_lifters)[i] for i in lifter_indices]
                step.append((weight, abs_lifter_indices))
                # Remove used lifters
                for idx in abs_lifter_indices:
                    available_lifters.remove(idx)
            else:
                remaining_boxes.append((weight, box_idx))
        
        steps.append(step)
        boxes = remaining_boxes
        
        if len(steps) > 7:  # Check if we exceed maximum steps
            return None
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

# Solve and print the result
print(solve_box_lifting())