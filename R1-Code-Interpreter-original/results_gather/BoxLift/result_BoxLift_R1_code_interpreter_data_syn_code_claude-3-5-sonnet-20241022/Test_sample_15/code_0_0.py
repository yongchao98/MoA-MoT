from itertools import combinations
import copy

def can_lift(lifters_combo, weight):
    return sum(lifters[i] for i in lifters_combo) >= weight

def find_solution():
    boxes = [50, 70, 83, 17, 21, 27, 53, 31]
    global lifters
    lifters = [42, 53, 74, 78]
    n_lifters = len(lifters)
    
    # Generate all possible lifter combinations
    all_lifter_combos = []
    for r in range(1, n_lifters + 1):
        all_lifter_combos.extend(combinations(range(n_lifters), r))
    
    def try_step(remaining_boxes, steps_so_far, max_steps=3):
        if not remaining_boxes:
            return steps_so_far
        if len(steps_so_far) >= max_steps:
            return None
            
        # Try all possible assignments for current step
        remaining_boxes = list(remaining_boxes)
        available_lifters = set(range(n_lifters))
        current_step = []
        
        # Sort boxes by weight (descending) for better heuristic
        box_indices = list(range(len(remaining_boxes)))
        box_indices.sort(key=lambda i: remaining_boxes[i], reverse=True)
        
        for box_idx in box_indices:
            if box_idx == -1:
                continue
            box_weight = remaining_boxes[box_idx]
            
            # Try all possible lifter combinations for this box
            for combo in all_lifter_combos:
                if not all(l in available_lifters for l in combo):
                    continue
                if can_lift(combo, box_weight):
                    # Try this assignment
                    current_step.append((box_weight, list(combo)))
                    for l in combo:
                        available_lifters.remove(l)
                    remaining_boxes[box_idx] = -1
                    break
        
        if not current_step:  # Couldn't make any assignments
            return None
            
        # Remove lifted boxes
        new_remaining = [b for b in remaining_boxes if b != -1]
        
        # Recurse
        new_steps = steps_so_far + [current_step]
        result = try_step(new_remaining, new_steps, max_steps)
        if result:
            return result
        return None

    solution = try_step(boxes, [])
    
    # Format output
    if solution:
        output = ""
        for i, step in enumerate(solution, 1):
            output += f"Step {i}: {step}\n"
        return output.strip()
    return "No solution found"

print(find_solution())