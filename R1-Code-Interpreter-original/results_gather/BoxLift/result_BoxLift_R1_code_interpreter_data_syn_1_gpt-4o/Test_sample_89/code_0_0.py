boxes = [227, 106, 388, 173, 63, 178, 261, 205, 323, 124, 386, 81, 71, 127, 134, 212, 150, 114, 41, 277, 123, 152, 47, 47]
lifters = [166, 106, 122, 181, 102, 119, 147]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        # Try to find a single lifter or a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                used_lifters[j] = True
                remaining_weight -= lifter
            if remaining_weight <= 0:
                break
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes[i] = None  # Mark this box as lifted
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")