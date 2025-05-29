from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Returns all possible combinations of lifters that can lift the weight
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                valid_combinations.append(list(combo))
    return valid_combinations

def solve_box_lifting():
    boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
    lifters = [78, 98, 51, 75, 43]
    
    # Sort boxes from heaviest to lightest
    boxes = sorted(boxes, reverse=True)
    
    # Store the solution
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 5:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in this step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find all possible lifter combinations for this box
            possible_lifters = [combo for combo in can_lift(lifters, box) 
                              if all(i in available_lifters for i in combo)]
            
            if possible_lifters:
                # Use the combination with minimum number of lifters
                best_combo = min(possible_lifters, key=len)
                step.append((box, best_combo))
                available_lifters -= set(best_combo)
                remaining_boxes.remove(box)
        
        if step:
            solution.append(step)
        else:
            break
    
    # Format and print solution
    if not remaining_boxes and len(solution) <= 5:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
        return True
    return False

solve_box_lifting()