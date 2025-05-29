from itertools import combinations
import copy

# Input data
boxes = [156, 123, 139, 191, 91, 193, 100, 47, 36, 169, 90, 117]
lifters = [52, 55, 63, 44]

def can_lift(weight, lifter_combo):
    return sum(lifters[i] for i in lifter_combo) >= weight

def find_solution():
    remaining_boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters
        for box_weight, box_idx in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                found = False
                for lifter_combo in combinations(available_lifters, r):
                    if can_lift(box_weight, lifter_combo):
                        step.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.remove((box_weight, box_idx))
                        found = True
                        break
                if found:
                    break
        
        if step:
            steps.append(step)
        else:
            return None
            
    return steps

# Get and format solution
solution = find_solution()
if solution:
    output = []
    for i, step in enumerate(solution, 1):
        output.append(f"Step {i}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found")