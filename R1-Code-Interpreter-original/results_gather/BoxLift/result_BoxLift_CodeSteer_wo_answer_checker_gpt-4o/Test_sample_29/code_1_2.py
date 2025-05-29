boxes = [100, 225, 68, 49, 127, 244, 299, 344, 330, 310, 299, 189, 124, 216, 350, 367, 85, 231, 165, 144, 323, 342, 368, 179]
lifters = [64, 93, 158, 137, 107, 95, 173]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    remaining_weight = box
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
        if remaining_weight <= 0:
            return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        steps.append(step.copy())
        return True
    if len(steps) >= 9:
        return False

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            for idx in lifter_indices:
                used[idx] = True
            step.append((box, lifter_indices))
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step, steps, used):
                return True
            step.pop()
            for idx in lifter_indices:
                used[idx] = False
    return False

steps = []
used = [False] * len(lifters)
backtrack(boxes, lifters, [], steps, used)

# Print the steps
output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"
output += ">>>"
print(output)