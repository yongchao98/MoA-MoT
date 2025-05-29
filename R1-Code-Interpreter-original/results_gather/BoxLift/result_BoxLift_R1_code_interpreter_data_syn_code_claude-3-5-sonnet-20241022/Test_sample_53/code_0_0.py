import copy
from collections import defaultdict

def can_lift(weight, lifters, capacity):
    if not lifters:
        return False
    return sum(capacity[i] for i in lifters) >= weight

def find_minimal_lifters(weight, available_lifters, capacity):
    if not available_lifters:
        return None
    
    # Try combinations of lifters from smallest to largest size
    for size in range(1, len(available_lifters) + 1):
        if size > 4:  # Limit combinations to prevent excessive computation
            break
        
        def try_combination(curr_lifters, start, remaining):
            if remaining == 0:
                if sum(capacity[i] for i in curr_lifters) >= weight:
                    return curr_lifters
                return None
            
            for i in range(start, len(available_lifters)):
                result = try_combination(curr_lifters + [available_lifters[i]], 
                                      i + 1, remaining - 1)
                if result:
                    return result
            return None
        
        result = try_combination([], 0, size)
        if result:
            return result
    return None

def solve_box_lifting(boxes, capacities):
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    steps = []
    remaining_boxes = set(range(len(boxes)))
    
    while remaining_boxes:
        current_step = []
        available_lifters = set(range(len(capacities)))
        
        # Try to assign as many boxes as possible in current step
        for box_idx, box_weight in boxes:
            if box_idx not in remaining_boxes or not available_lifters:
                continue
                
            # Find minimal set of lifters that can lift this box
            lifters = find_minimal_lifters(box_weight, list(available_lifters), capacities)
            
            if lifters:
                current_step.append((box_weight, lifters))
                remaining_boxes.remove(box_idx)
                available_lifters -= set(lifters)
        
        if current_step:
            steps.append(current_step)
        else:
            break
            
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

# Input data
boxes = [108, 206, 326, 343, 45, 374, 227, 366, 110, 136, 343, 298, 201, 168, 191, 243, 365, 375, 294, 317]
lifters = [147, 80, 144, 133, 69, 200]

solve_box_lifting(boxes, lifters)