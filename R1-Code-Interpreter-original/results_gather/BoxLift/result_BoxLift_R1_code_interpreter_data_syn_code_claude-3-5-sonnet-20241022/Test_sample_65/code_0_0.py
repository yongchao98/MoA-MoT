from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Try all possible combinations of lifters to see if they can lift the weight
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [70, 51, 67, 51, 34, 82, 20, 93]
    lifter_capacities = [64, 51, 79, 80]
    boxes_with_indices = list(enumerate(boxes))
    solution = []
    
    while boxes_with_indices:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort boxes by weight (descending) to handle heaviest boxes first
        boxes_with_indices.sort(key=lambda x: x[1], reverse=True)
        
        # Try to assign boxes to lifters
        remaining_boxes = []
        for box_idx, weight in boxes_with_indices:
            if not available_lifters:
                remaining_boxes.append((box_idx, weight))
                continue
                
            # Find lifters that can handle this box
            possible_lifters = []
            for combo_size in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, combo_size):
                    if sum(lifter_capacities[i] for i in combo) >= weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step.append((weight, list(possible_lifters)))
                available_lifters -= set(possible_lifters)
            else:
                remaining_boxes.append((box_idx, weight))
        
        boxes_with_indices = remaining_boxes
        if step:
            solution.append(step)
            
        if len(solution) > 3:
            return None
    
    return solution

# Get and format solution
solution = solve_box_lifting()
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(output.strip())
else:
    print("No solution found within 3 steps")