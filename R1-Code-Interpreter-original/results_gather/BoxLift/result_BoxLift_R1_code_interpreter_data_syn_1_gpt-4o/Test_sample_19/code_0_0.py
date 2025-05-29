boxes = [55, 77, 274, 297, 70, 312, 108, 84, 340, 273, 157, 386, 95, 223, 311, 306, 287, 91, 220, 296, 320, 237, 97, 73]
lifters = [167, 138, 133, 87, 142, 167, 140]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        # Try to find a combination of lifters to lift the box
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
    boxes = [b for b in boxes if b is not None]  # Remove lifted boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")