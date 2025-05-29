boxes = [75, 298, 70, 200, 278, 35, 278, 56, 270, 212, 161, 98, 244, 292, 177, 145]
lifters = [105, 102, 125, 93, 123]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters that can lift the box
    lifter_indices = []
    total_capacity = 0
    for i, capacity in enumerate(lifters):
        if not used[i] and total_capacity + capacity <= box:
            total_capacity += capacity
            lifter_indices.append(i)
        if total_capacity >= box:
            return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 7:
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
    print("No solution found within 7 steps.")