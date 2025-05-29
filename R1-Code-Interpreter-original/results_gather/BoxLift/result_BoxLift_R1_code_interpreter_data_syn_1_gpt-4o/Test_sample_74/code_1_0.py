from itertools import combinations

boxes = [166, 144, 53, 213, 51, 156, 197, 311, 177, 358, 172, 134, 179, 145, 91, 188, 352, 294, 292, 88, 97, 394, 123, 294]
lifters = [88, 185, 145, 195, 147, 145]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 8:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        if box is None:
            continue
        # Try to find a single lifter or a combination of lifters to lift the box
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
    print(f"Step {i+1}: {step}")