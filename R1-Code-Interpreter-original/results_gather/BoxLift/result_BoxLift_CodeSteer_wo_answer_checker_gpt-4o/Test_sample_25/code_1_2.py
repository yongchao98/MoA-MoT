def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    
    for i, capacity in enumerate(lifters):
        if not used[i] and remaining_weight > 0:
            if capacity <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= capacity
                used[i] = True
    
    if remaining_weight <= 0:
        return lifter_indices
    else:
        for i in lifter_indices:
            used[i] = False
        return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    
    if step >= 4:
        return False
    
    current_step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if backtrack(remaining_boxes, lifters, step + 1, steps, [False] * len(lifters)):
            return True
        steps.pop()
    
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    used = [False] * len(lifters)
    
    if backtrack(boxes, lifters, 0, steps, used):
        return steps
    else:
        return "Not possible to lift all boxes in 4 steps or less."

boxes = [17, 69, 16, 83, 95, 91, 68, 25]
lifters = [42, 68, 50, 48]

steps = assign_lifters_to_boxes(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)