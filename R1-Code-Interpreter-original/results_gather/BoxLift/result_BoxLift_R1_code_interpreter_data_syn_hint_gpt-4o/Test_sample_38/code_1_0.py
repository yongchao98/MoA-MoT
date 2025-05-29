boxes = [356, 292, 240, 155, 142, 223, 178, 231, 254, 315, 236, 328, 203, 397, 44, 200, 80, 131, 175, 216]
lifters = [155, 116, 53, 98, 96, 156, 167]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
max_steps = 7

def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                return lifter_indices
    return None

for _ in range(max_steps):
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []

    for box in boxes:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                used_lifters[index] = True
        else:
            remaining_boxes.append(box)

    steps.append(step)
    boxes = remaining_boxes

    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")