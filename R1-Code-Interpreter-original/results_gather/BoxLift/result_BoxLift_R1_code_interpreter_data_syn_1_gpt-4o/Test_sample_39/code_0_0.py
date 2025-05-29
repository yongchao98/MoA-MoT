boxes = [256, 53, 286, 205, 261, 209, 242, 143, 354, 55, 72, 161, 268, 305, 43, 246, 146, 331, 84, 128, 88, 358, 71, 93]
lifters = [189, 58, 58, 79, 123, 192, 116]

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
                used_lifters[j] = True
                remaining_weight -= lifter
            if remaining_weight <= 0:
                break
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
    # Remove lifted boxes
    for box, _ in step:
        boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")