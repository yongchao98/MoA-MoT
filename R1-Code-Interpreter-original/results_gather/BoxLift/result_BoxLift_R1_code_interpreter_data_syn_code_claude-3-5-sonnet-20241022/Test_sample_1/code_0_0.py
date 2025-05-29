from itertools import combinations
import copy

def can_lift_box(lifters_indices, lifters_capacity, box_weight):
    return sum(lifters_capacity[i] for i in lifters_indices) >= box_weight

def find_solution():
    boxes = [71, 58, 105, 85, 123, 125, 84, 105, 141, 73, 102, 46]
    lifters = [77, 45, 120, 109]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes and len(steps) < 5:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        current_boxes = remaining_boxes.copy()
        for box_idx, box_weight in current_boxes:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            assigned = False
            for r in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifter_combo in combinations(available_lifters, r):
                    if can_lift_box(lifter_combo, lifters, box_weight):
                        step.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.remove((box_idx, box_weight))
                        assigned = True
                        break
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format and print solution
    if not remaining_boxes and len(steps) <= 5:
        result = ""
        for i, step in enumerate(steps, 1):
            result += f"Step {i}: {step}\n"
        return result.strip()
    return "No solution found within 5 steps"

print(find_solution())