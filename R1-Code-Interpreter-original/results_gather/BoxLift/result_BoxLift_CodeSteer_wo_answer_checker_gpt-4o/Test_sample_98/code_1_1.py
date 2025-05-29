def can_lift(box, lifters, used):
    # Try to lift the box with available lifters
    lifter_indices = []
    remaining_weight = box
    
    for i, capacity in enumerate(lifters):
        if not used[i] and capacity >= remaining_weight:
            lifter_indices.append(i)
            used[i] = True
            return lifter_indices
    
    for i, capacity in enumerate(lifters):
        if not used[i] and capacity <= remaining_weight:
            lifter_indices.append(i)
            used[i] = True
            remaining_weight -= capacity
            if remaining_weight <= 0:
                return lifter_indices
    
    # If we can't lift the box with available lifters
    for index in lifter_indices:
        used[index] = False
    return None

def backtrack(boxes, lifters, step, steps, max_steps):
    if not boxes:
        return True
    
    if step >= max_steps:
        return False
    
    used = [False] * len(lifters)
    current_step = []
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if backtrack(boxes, lifters, step + 1, steps, max_steps):
            return True
        steps.pop()
        boxes.extend([box for box, _ in current_step])
    
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    steps = []
    max_steps = 6
    
    if backtrack(boxes, lifters, 0, steps, max_steps):
        return steps
    else:
        return None

boxes = [55, 244, 173, 293, 90, 126, 340, 250, 66, 143, 103, 244, 76, 166, 130, 216, 54, 196, 245, 307]
lifters = [142, 178, 196, 52, 101, 144, 50]

steps = assign_lifters_to_boxes(boxes, lifters)

if steps:
    output = ""
    for i, step in enumerate(steps):
        output += f"Step {i + 1}: {step}\n"
    print(f"<<<{output}>>>")
else:
    print("<<<No valid solution found within 6 steps>>>")