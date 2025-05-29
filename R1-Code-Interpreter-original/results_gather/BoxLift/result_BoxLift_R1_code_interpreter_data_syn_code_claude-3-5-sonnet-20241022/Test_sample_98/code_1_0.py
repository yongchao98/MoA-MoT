from itertools import combinations

def can_lift(lifters_combo, weight, lifter_capacities):
    return sum(lifter_capacities[i] for i in lifters_combo) >= weight

def find_solution(boxes, lifter_capacities):
    # Sort boxes in descending order and create index mapping
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    steps = []
    
    while boxes_with_index:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Process boxes in current step
        boxes_to_remove = []
        for box_idx, box_weight in boxes_with_index:
            if not available_lifters:  # No more lifters available
                break
                
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                lifter_combos = combinations(available_lifters, r)
                assigned = False
                for lifters_combo in lifter_combos:
                    if can_lift(lifters_combo, box_weight, lifter_capacities):
                        step.append((box_weight, list(lifters_combo)))
                        available_lifters -= set(lifters_combo)
                        boxes_to_remove.append((box_idx, box_weight))
                        assigned = True
                        break
                if assigned:
                    break
        
        # Remove processed boxes
        for box in boxes_to_remove:
            boxes_with_index.remove(box)
            
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