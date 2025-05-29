import itertools

def can_lift(lifters, weight):
    # Check all possible combinations of lifters to see if they can lift the weight
    for r in range(1, len(lifters) + 1):
        for combo in itertools.combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [45, 139, 291, 189, 171, 265, 139, 312, 289, 394, 183, 229, 321, 250, 226, 335, 59, 127, 298, 339, 256, 161, 288, 179]
    lifters = [52, 99, 86, 56, 92, 134]
    
    # Sort boxes from heaviest to lightest
    boxes_with_index = list(enumerate(boxes))
    boxes_with_index.sort(key=lambda x: x[1], reverse=True)
    
    remaining_boxes = boxes_with_index.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to lifters
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try to find lifters that can lift this box
            possible_lifters = []
            for r in range(1, len(available_lifters) + 1):
                for combo in itertools.combinations(available_lifters, r):
                    if sum(lifters[i] for i in combo) >= box_weight:
                        possible_lifters = list(combo)
                        break
                if possible_lifters:
                    break
            
            if possible_lifters:
                step.append((box_weight, possible_lifters))
                available_lifters -= set(possible_lifters)
                remaining_boxes.remove((box_idx, box_weight))
        
        steps.append(step)
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()