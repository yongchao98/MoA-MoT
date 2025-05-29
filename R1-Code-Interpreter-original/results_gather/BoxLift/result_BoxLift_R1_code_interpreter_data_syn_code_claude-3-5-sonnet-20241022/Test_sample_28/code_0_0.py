from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [71, 72, 167, 64, 339, 45, 243, 53, 52, 280, 164, 322, 362, 75, 158, 142, 150, 353, 182, 282]
    lifters = [139, 55, 172, 82, 129, 59, 130]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    
    solution = []
    remaining_boxes = boxes[:]
    
    while remaining_boxes and len(solution) < 7:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        for box_weight, box_idx in remaining_boxes[:]:
            # Find lifters who can lift this box
            possible_lifters = []
            for r in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, r):
                    if sum(lifters[i] for i in combo) >= box_weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step.append((box_weight, possible_lifters))
                available_lifters -= set(possible_lifters)
                remaining_boxes.remove((box_weight, box_idx))
            
            if not available_lifters:
                break
        
        if step:
            solution.append(step)
        else:
            break
    
    # Print solution
    if not remaining_boxes and len(solution) <= 7:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
    else:
        print("No solution found within 7 steps!")

solve_box_lifting()