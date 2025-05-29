def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    
    for i, capacity in enumerate(lifters):
        if not used[i] and remaining_weight > 0:
            lifter_indices.append(i)
            remaining_weight -= capacity
            used[i] = True
            if remaining_weight <= 0:
                break
    
    if remaining_weight > 0:
        for i in lifter_indices:
            used[i] = False
        return False, []
    
    return True, lifter_indices

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    
    if step >= 7:
        return False
    
    current_step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        can_lift_box, lifter_indices = can_lift(box, lifters, used)
        if can_lift_box:
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
        return []

boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

steps = assign_lifters_to_boxes(boxes, lifters)

output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output}>>>")