from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Try combinations of lifters from 1 to len(lifters)
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [166, 194, 27, 144, 198, 172, 81, 117, 154, 183, 27, 184]
    lifter_capacities = [97, 90, 80, 111, 76]
    
    # Sort boxes from heaviest to lightest
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    steps = []
    remaining_boxes = boxes_with_index.copy()
    
    while remaining_boxes and len(steps) < 5:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign as many boxes as possible in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find lifters who can handle this box
            possible_lifters = []
            for i in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, i):
                    if sum(lifter_capacities[j] for j in combo) >= box_weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step.append((box_weight, possible_lifters))
                available_lifters -= set(possible_lifters)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    if not remaining_boxes and len(steps) <= 5:
        print(f"<<<{result.strip()}>>>")
    else:
        print("No solution found within 5 steps")

solve_box_lifting()