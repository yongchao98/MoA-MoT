boxes = [12, 90, 67, 10, 10, 94, 60, 46]
lifters = [56, 77, 72]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    lifter_indices = []
    total_capacity = 0
    for i, capacity in enumerate(lifters):
        if not used[i] and total_capacity < box:
            lifter_indices.append(i)
            total_capacity += capacity
            used[i] = True
        if total_capacity >= box:
            return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        steps.append(step.copy())
        return True

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step, steps, used):
                return True
            step.pop()
            for index in lifter_indices:
                used[index] = False

    return False

steps = []
used = [False] * len(lifters)
backtrack(boxes, lifters, [], steps, used)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")