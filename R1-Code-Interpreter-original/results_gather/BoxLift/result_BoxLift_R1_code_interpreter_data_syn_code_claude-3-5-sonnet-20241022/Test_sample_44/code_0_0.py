import copy
from collections import defaultdict

def can_lift(lifters, weight):
    if not lifters:
        return False
    return sum(lifters) >= weight

def find_lifter_combination(available_lifters, weight, max_lifters=6):
    if not available_lifters:
        return None
    
    # Try combinations from smallest to largest size
    available_lifters = sorted(available_lifters, reverse=True)
    n = len(available_lifters)
    
    # Try with minimum number of lifters first
    for size in range(1, min(max_lifters + 1, n + 1)):
        # Try all combinations of given size
        def try_combinations(curr_sum, curr_indices, start, size):
            if size == 0:
                if curr_sum >= weight:
                    return curr_indices
                return None
            
            for i in range(start, n):
                result = try_combinations(curr_sum + available_lifters[i], 
                                       curr_indices + [i], i + 1, size - 1)
                if result:
                    return result
            return None
        
        result = try_combinations(0, [], 0, size)
        if result:
            return result
    return None

def solve_box_lifting():
    boxes = [116, 246, 296, 369, 78, 275, 77, 383, 71, 155, 352, 94, 63, 168, 
             350, 79, 59, 252, 88, 278, 188, 383, 240, 308]
    lifter_capacities = [134, 132, 150, 137, 62, 132]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    steps = []
    remaining_boxes = boxes[:]
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifter_capacities)))
        current_boxes = remaining_boxes[:]
        
        while current_boxes and available_lifters:
            box_weight, box_index = current_boxes[0]
            
            # Find lifters combination for current box
            lifter_indices = find_lifter_combination(
                [lifter_capacities[i] for i in available_lifters],
                box_weight
            )
            
            if lifter_indices is not None:
                # Convert relative indices to absolute indices
                abs_lifter_indices = [available_lifters[i] for i in lifter_indices]
                step.append((box_weight, abs_lifter_indices))
                
                # Remove used lifters
                for idx in sorted(lifter_indices, reverse=True):
                    available_lifters.pop(idx)
                
                # Remove assigned box
                remaining_boxes.remove((box_weight, box_index))
            
            current_boxes.pop(0)
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())