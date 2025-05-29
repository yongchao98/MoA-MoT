boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        if not boxes:
            break
        # Try to find a combination of lifters to lift the current box
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
            boxes[i] = None  # Mark this box as lifted
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)
    if len(steps) > 7:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")