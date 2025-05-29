from itertools import combinations
import copy

def can_lift(weight, lifters, capacities):
    return sum(capacities[i] for i in lifters) >= weight

def find_valid_lifter_combinations(weight, available_lifters, capacities):
    valid_combinations = []
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if can_lift(weight, combo, capacities):
                valid_combinations.append(list(combo))
    return valid_combinations

def solve_box_lifting():
    boxes = [86, 183, 78, 85, 112, 170, 130, 121, 208, 226, 39, 176, 68, 256, 56, 34]
    capacities = [116, 81, 149, 77, 138]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    step = 0
    
    while remaining_boxes and step < 6:
        step += 1
        step_solution = []
        available_lifters = set(range(len(capacities)))
        
        # Try to assign boxes to lifters
        boxes_to_remove = []
        for box in remaining_boxes:
            if not available_lifters:
                break
                
            valid_combinations = find_valid_lifter_combinations(box, list(available_lifters), capacities)
            if valid_combinations:
                # Use the combination with minimum number of lifters
                best_combo = min(valid_combinations, key=len)
                step_solution.append((box, best_combo))
                available_lifters -= set(best_combo)
                boxes_to_remove.append(box)
        
        # Remove assigned boxes
        for box in boxes_to_remove:
            remaining_boxes.remove(box)
            
        solution.append(step_solution)
    
    # Format the solution
    formatted_solution = []
    for i, step in enumerate(solution, 1):
        formatted_solution.append(f"Step {i}: {step}")
    
    return "\n".join(formatted_solution) if not remaining_boxes else "No solution found"

print(solve_box_lifting())