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
    boxes = [97, 164, 187, 145, 255, 257, 240, 55, 261, 185, 252, 257, 206, 298, 121, 202]
    lifters = [131, 106, 144, 121, 109, 52]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    steps = []
    
    while boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        i = 0
        while i < len(boxes):
            box_weight, box_index = boxes[i]
            
            # Try different combinations of available lifters
            best_lifters = None
            for r in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, r):
                    if sum(lifters[j] for j in combo) >= box_weight:
                        best_lifters = list(combo)
                        break
                if best_lifters:
                    break
            
            if best_lifters:
                step.append((box_weight, best_lifters))
                available_lifters -= set(best_lifters)
                boxes.pop(i)
            else:
                i += 1
                
            if not available_lifters:
                break
        
        steps.append(step)
        
        if len(steps) > 7:
            return None
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

solution = solve_box_lifting()
print(solution)