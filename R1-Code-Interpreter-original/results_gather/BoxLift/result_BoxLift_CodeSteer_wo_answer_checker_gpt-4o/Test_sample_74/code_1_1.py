boxes = [166, 144, 53, 213, 51, 156, 197, 311, 177, 358, 172, 134, 179, 145, 91, 188, 352, 294, 292, 88, 97, 394, 123, 294]
lifters = [88, 185, 145, 195, 147, 145]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    # Try to find a combination of lifters to lift the box
    lifter_indices = []
    remaining_weight = box
    
    for i, lifter in enumerate(lifters):
        if not used_lifters[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            used_lifters[i] = True
        
        if remaining_weight <= 0:
            return lifter_indices
    
    # If we can't lift the box with available lifters, backtrack
    for i in lifter_indices:
        used_lifters[i] = False
    
    return None

def backtrack(boxes, lifters, step, steps, used_boxes, max_steps):
    if step >= max_steps:
        return False
    
    if all(used_boxes):
        return True
    
    used_lifters = [False] * len(lifters)
    current_step = []
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        lifter_indices = can_lift(box, lifters, used_lifters)
        
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            used_boxes[i] = True
    
    if current_step:
        steps.append(current_step)
        if backtrack(boxes, lifters, step + 1, steps, used_boxes, max_steps):
            return True
        steps.pop()
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            used_boxes[i] = False
    
    return False

steps = []
used_boxes = [False] * len(boxes)
max_steps = 8

if backtrack(boxes, lifters, 0, steps, used_boxes, max_steps):
    output = []
    for idx, step in enumerate(steps):
        step_str = f"Step {idx + 1}: " + str(step)
        output.append(step_str)
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within the given step limit.")