from itertools import combinations

def can_lift(box_weight, lifters_combo):
    return sum(lifters[i] for i in lifters_combo) >= box_weight

def find_solution():
    boxes = [92, 26, 84, 55, 72, 48, 33, 19]
    global lifters
    lifters = [64, 48, 63]
    boxes_indices = list(range(len(boxes)))
    remaining_boxes = list(zip(boxes, boxes_indices))
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Sort boxes by weight (descending)
        remaining_boxes.sort(reverse=True)
        
        # Try to assign boxes to lifters
        for box_weight, box_idx in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                found = False
                for lifters_combo in combinations(available_lifters, r):
                    if can_lift(box_weight, lifters_combo):
                        step.append((box_weight, list(lifters_combo)))
                        available_lifters -= set(lifters_combo)
                        remaining_boxes.remove((box_weight, box_idx))
                        found = True
                        break
                if found:
                    break
        
        steps.append(step)
        if len(steps) > 4:  # Check if we exceed 4 steps
            return None
    
    return steps

# Find and format solution
solution = find_solution()
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(output.strip())
else:
    print("No solution found within 4 steps")