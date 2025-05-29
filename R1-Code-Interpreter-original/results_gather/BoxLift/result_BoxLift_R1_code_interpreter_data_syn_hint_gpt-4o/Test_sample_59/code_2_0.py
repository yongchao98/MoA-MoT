boxes = [299, 338, 270, 207, 76, 195, 101, 89, 277, 236, 240, 342, 47, 207, 220, 349, 185, 124, 309, 112, 266, 278, 162, 119]
lifters = [198, 55, 76, 137, 163, 117, 146]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 8:
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
            boxes[i] = None  # Mark this box as lifted
    boxes = [b for b in boxes if b is not None]  # Remove lifted boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")