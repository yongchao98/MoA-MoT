boxes = [110, 230, 255, 47, 133, 280, 271, 275, 155, 80, 169, 89, 299, 241, 187, 234]
lifters = [160, 110, 114, 55, 54, 60]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters to lift the box
    lifter_indices = []
    remaining_weight = box
    
    for j, lifter in enumerate(lifters):
        if not used[j] and lifter <= remaining_weight:
            lifter_indices.append(j)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                return lifter_indices
    
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    
    if step >= 7:
        return False
    
    current_step = []
    for i, box in enumerate(boxes[:]):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            for idx in lifter_indices:
                used[idx] = True
            boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters)):
            return True
        steps.pop()
        boxes.extend([b[0] for b in current_step])
    
    return False

steps = []
used = [False] * len(lifters)
backtrack(boxes, lifters, 0, steps, used)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")