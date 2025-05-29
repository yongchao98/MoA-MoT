from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Returns list of lifter combinations that can lift the weight
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                valid_combinations.append(list(combo))
    return valid_combinations

def solve_box_lifting():
    boxes = [187, 170, 257, 238, 229, 202, 239, 125, 197, 243, 110, 50, 48, 168, 42, 292]
    lifters = [115, 55, 142, 79, 63]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find valid lifter combinations for this box
            valid_combos = [c for c in can_lift(lifters, box) 
                          if all(i in available_lifters for i in c)]
            
            if valid_combos:
                # Use the combination with minimum number of lifters
                best_combo = min(valid_combos, key=len)
                step.append((box, best_combo))
                available_lifters -= set(best_combo)
                remaining_boxes.remove(box)
        
        if step:
            solution.append(step)
        else:
            # If we can't make progress, something is wrong
            return None
            
        if len(solution) > 10:
            return None
    
    return solution

# Solve and format output
solution = solve_box_lifting()
if solution:
    output = []
    for i, step in enumerate(solution, 1):
        step_str = f"Step {i}: {step}"
        output.append(step_str)
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("<<<No solution found within 10 steps>>>")