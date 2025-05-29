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
    boxes = [24, 78, 95, 23, 57, 36, 40, 26]
    lifter_capacities = [56, 51, 66, 65]
    max_steps = 3
    
    def try_assignment(remaining_boxes, steps_so_far):
        if not remaining_boxes:
            return steps_so_far
        if len(steps_so_far) >= max_steps:
            return None
            
        # Try all possible assignments for current step
        current_step = []
        available_lifters = set(range(len(lifter_capacities)))
        remaining_boxes_copy = remaining_boxes.copy()
        
        while remaining_boxes_copy and available_lifters:
            assigned = False
            for i, box in enumerate(remaining_boxes_copy):
                # Try to lift current box with available lifters
                temp_lifters = [lifter_capacities[j] for j in available_lifters]
                lifter_indices = can_lift([lifter_capacities[j] for j in available_lifters], box)
                
                if lifter_indices is not None:
                    # Convert relative indices to absolute indices
                    abs_indices = [list(available_lifters)[j] for j in lifter_indices]
                    current_step.append((box, abs_indices))
                    # Remove used lifters
                    for idx in lifter_indices:
                        available_lifters.remove(list(available_lifters)[idx])
                    remaining_boxes_copy.pop(i)
                    assigned = True
                    break
            
            if not assigned:
                break
        
        if current_step:
            # Remove assigned boxes from remaining boxes
            new_remaining = [b for b in remaining_boxes if not any(b == step[0] for step in current_step)]
            result = try_assignment(new_remaining, steps_so_far + [current_step])
            if result is not None:
                return result
        
        return None

    solution = try_assignment(boxes, [])
    
    if solution:
        for i, step in enumerate(solution, 1):
            print(f"Step {i}: {step}")
    else:
        print("No solution found within 3 steps")

solve_box_lifting()