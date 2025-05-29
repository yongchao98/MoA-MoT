import itertools

def can_lift(box_weight, lifters_combination, lifter_capacities):
    return sum(lifter_capacities[i] for i in lifters_combination) >= box_weight

def find_solution():
    boxes = [169, 198, 165, 154, 71, 159, 55, 205, 299, 170, 122, 160, 43, 259, 246, 172]
    lifter_capacities = [80, 51, 45, 141, 83, 152]
    
    # Sort boxes in descending order
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            assigned = False
            for size in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifters in itertools.combinations(available_lifters, size):
                    if can_lift(box_weight, lifters, lifter_capacities):
                        step.append((box_weight, list(lifters)))
                        available_lifters -= set(lifters)
                        remaining_boxes.remove((box_idx, box_weight))
                        assigned = True
                        break
        
        if step:
            solution.append(step)
        
        if len(solution) > 7:
            return None
    
    return solution

# Find and format solution
solution = find_solution()
if solution:
    formatted_solution = []
    for i, step in enumerate(solution, 1):
        formatted_solution.append(f"Step {i}: {step}")
    print("<<<" + "\n".join(formatted_solution) + ">>>")
else:
    print("<<<No solution found within 7 steps>>>")