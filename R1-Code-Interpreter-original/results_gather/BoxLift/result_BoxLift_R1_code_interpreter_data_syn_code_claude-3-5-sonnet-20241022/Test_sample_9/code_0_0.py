from itertools import combinations

def can_lift(box_weight, lifters_combo):
    return sum(lifters[i] for i in lifters_combo) >= box_weight

def solve_box_lifting():
    boxes = [355, 40, 260, 256, 250, 107, 301, 298, 72, 316, 50, 319, 243, 51, 393, 268, 141, 178, 270, 158, 354, 336, 55, 84]
    global lifters
    lifters = [100, 90, 66, 141, 184, 78, 50]
    
    # Sort boxes from heaviest to lightest
    boxes = sorted(boxes, reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in this step
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Try different combinations of available lifters
            found_lifters = None
            for size in range(1, len(available_lifters) + 1):
                if found_lifters:
                    break
                for lifter_combo in combinations(available_lifters, size):
                    if can_lift(box, lifter_combo):
                        found_lifters = lifter_combo
                        break
            
            if found_lifters:
                step.append((box, list(found_lifters)))
                available_lifters -= set(found_lifters)
                remaining_boxes.remove(box)
        
        steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(f"<<{result}>>")

solve_box_lifting()