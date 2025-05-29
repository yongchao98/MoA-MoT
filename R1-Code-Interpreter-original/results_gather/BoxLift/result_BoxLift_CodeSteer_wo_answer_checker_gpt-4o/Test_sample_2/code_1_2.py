boxes = [66, 299, 90, 133, 52, 245, 57, 62, 71, 217, 117, 193, 204, 84, 224, 51]
lifters = [123, 155, 65, 92, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a single lifter or a combination of lifters to lift the box
    for i in range(len(lifters)):
        if used[i]:
            continue
        if lifters[i] >= box:
            return [i]
    for i in range(len(lifters)):
        if used[i]:
            continue
        for j in range(i + 1, len(lifters)):
            if used[j]:
                continue
            if lifters[i] + lifters[j] >= box:
                return [i, j]
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 6:
        return False

    current_step = []
    remaining_boxes = boxes[:]
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            for index in lifter_indices:
                used[index] = True
            remaining_boxes.remove(box)

    if current_step:
        steps.append(current_step)
        if backtrack(remaining_boxes, lifters, step + 1, steps, [False] * len(lifters)):
            return True
        steps.pop()

    return False

steps = []
used = [False] * len(lifters)
if backtrack(boxes, lifters, 0, steps, used):
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 6 steps.")