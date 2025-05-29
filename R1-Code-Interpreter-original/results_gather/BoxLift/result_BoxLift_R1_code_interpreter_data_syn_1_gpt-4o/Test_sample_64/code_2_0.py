boxes = [45, 139, 291, 189, 171, 265, 139, 312, 289, 394, 183, 229, 321, 250, 226, 335, 59, 127, 298, 339, 256, 161, 288, 179]
lifters = [52, 99, 86, 56, 92, 134]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        lifter_indices = []
        remaining_weight = box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                remaining_weight -= lifter
                used_lifters[j] = True
            if remaining_weight <= 0:
                break
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes[i] = None  # Mark the box as lifted
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")