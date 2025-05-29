from itertools import combinations
import copy

def can_lift(weight, lifters_capacity):
    return sum(lifters_capacity) >= weight

def find_solution():
    boxes = [66, 299, 90, 133, 52, 245, 57, 62, 71, 217, 117, 193, 204, 84, 224, 51]
    lifters = [123, 155, 65, 92, 92]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifters)))
        current_boxes = remaining_boxes.copy()
        
        while available_lifters and current_boxes:
            box_idx, box_weight = current_boxes[0]
            
            # Try combinations of available lifters
            best_lifters = None
            for r in range(1, len(available_lifters) + 1):
                for lifter_combo in combinations(available_lifters, r):
                    lifter_capacity = sum(lifters[i] for i in lifter_combo)
                    if lifter_capacity >= box_weight:
                        best_lifters = list(lifter_combo)
                        break
                if best_lifters:
                    break
            
            if best_lifters:
                step.append((box_weight, best_lifters))
                for lifter in best_lifters:
                    available_lifters.remove(lifter)
                current_boxes.pop(0)
                remaining_boxes.remove((box_idx, box_weight))
            else:
                current_boxes.pop(0)
                
        solution.append(step)
        
    # Format output
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    
    print(output.strip())

find_solution()