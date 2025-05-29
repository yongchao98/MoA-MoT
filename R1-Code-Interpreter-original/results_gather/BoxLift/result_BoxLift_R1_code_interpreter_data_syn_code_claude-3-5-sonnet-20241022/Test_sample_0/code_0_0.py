from itertools import combinations

def can_lift_box(box_weight, lifters_combo, lifter_capacities):
    return sum(lifter_capacities[i] for i in lifters_combo) >= box_weight

def find_solution():
    boxes = [18, 93, 23, 70, 72, 22, 96, 96]
    lifter_capacities = [60, 46, 61]
    num_lifters = len(lifter_capacities)
    
    # Keep track of remaining boxes
    remaining_boxes = list(enumerate(boxes))  # (index, weight)
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(num_lifters))
        
        # Try to assign boxes to lifters
        remaining_boxes.sort(key=lambda x: x[1], reverse=True)  # Sort by weight
        
        # For each remaining box
        i = 0
        while i < len(remaining_boxes):
            box_idx, box_weight = remaining_boxes[i]
            assigned = False
            
            # Try different combinations of available lifters
            for combo_size in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                    
                for lifter_combo in combinations(available_lifters, combo_size):
                    if can_lift_box(box_weight, lifter_combo, lifter_capacities):
                        step.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.pop(i)
                        assigned = True
                        break
            
            if not assigned:
                i += 1
                
        steps.append(step)
        
        if len(steps) > 5:  # If solution exceeds 5 steps
            return None
            
    return steps

# Find and format solution
solution = find_solution()
if solution:
    formatted_solution = []
    for i, step in enumerate(solution, 1):
        formatted_solution.append(f"Step {i}: {step}")
    print("\n".join(formatted_solution))
else:
    print("No solution found within 5 steps")