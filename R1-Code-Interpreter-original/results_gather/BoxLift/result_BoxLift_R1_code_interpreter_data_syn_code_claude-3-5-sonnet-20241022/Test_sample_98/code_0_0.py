from itertools import combinations

def can_lift(lifters_combo, weight, lifter_capacities):
    return sum(lifter_capacities[i] for i in lifters_combo) >= weight

def find_solution(boxes, lifter_capacities):
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    boxes_with_index = list(enumerate(boxes))
    steps = []
    
    while boxes_with_index:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters in this step
        remaining_boxes = boxes_with_index.copy()
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            box_idx, box_weight = remaining_boxes[i]
            
            # Try different combinations of available lifters
            assigned = False
            for r in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifters_combo in combinations(available_lifters, r):
                    if can_lift(lifters_combo, box_weight, lifter_capacities):
                        step.append((box_weight, list(lifters_combo)))
                        available_lifters -= set(lifters_combo)
                        boxes_with_index.remove((box_idx, box_weight))
                        assigned = True
                        break
            
            if not assigned:
                i += 1
                
        steps.append(step)
        
        if len(steps) > 6:  # If we exceed 6 steps, return None
            return None
            
    return steps

# Input data
boxes = [55, 244, 173, 293, 90, 126, 340, 250, 66, 143, 103, 244, 76, 166, 130, 216, 54, 196, 245, 307]
lifter_capacities = [142, 178, 196, 52, 101, 144, 50]

# Find solution
solution = find_solution(boxes, lifter_capacities)

# Format and print solution
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(output.strip())
else:
    print("No solution found within 6 steps")