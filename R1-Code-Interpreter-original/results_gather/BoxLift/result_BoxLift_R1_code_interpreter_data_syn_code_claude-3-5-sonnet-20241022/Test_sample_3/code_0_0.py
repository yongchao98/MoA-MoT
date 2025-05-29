from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for comb in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in comb) >= weight:
                return list(comb)
    return None

def solve_box_lifting():
    boxes = [159, 197, 211, 302, 51, 251, 299, 229, 384, 196, 377, 283, 127, 312, 228, 280, 153, 199, 169, 151]
    lifters = [52, 80, 175, 75, 136, 143]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    steps = []
    
    while boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        remaining_boxes = []
        for box_weight, box_idx in boxes:
            if not available_lifters:
                remaining_boxes.append((box_weight, box_idx))
                continue
                
            # Find lifters that can lift this box
            lifters_subset = [lifters[i] for i in available_lifters]
            lifter_indices = can_lift([lifters[i] for i in available_lifters], box_weight)
            
            if lifter_indices is not None:
                # Convert relative indices to absolute indices
                abs_indices = [list(available_lifters)[i] for i in lifter_indices]
                step.append((box_weight, abs_indices))
                available_lifters -= set(abs_indices)
            else:
                remaining_boxes.append((box_weight, box_idx))
        
        steps.append(step)
        boxes = remaining_boxes
        
        if len(steps) > 10:
            return "No solution found within 10 steps"
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())