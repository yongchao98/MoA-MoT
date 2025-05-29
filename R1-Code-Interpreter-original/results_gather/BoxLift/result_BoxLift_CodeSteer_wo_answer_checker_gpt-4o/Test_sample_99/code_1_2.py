boxes = [289, 375, 107, 145, 257, 48, 141, 83, 136, 368, 59, 133, 186, 266, 353, 73, 66, 210, 247, 79, 342, 318, 337, 162]
lifters = [75, 189, 162, 181, 118, 194, 137]

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
            used[j] = True
        if remaining_weight <= 0:
            return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 7:
        return False

    current_step = []
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            boxes.remove(box)

    if current_step:
        steps.append(current_step)
        if backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters)):
            return True
        steps.pop()
        boxes.extend([box for box, _ in current_step])

    return False

steps = []
used = [False] * len(lifters)
if backtrack(boxes, lifters, 0, steps, used):
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 7 steps.")