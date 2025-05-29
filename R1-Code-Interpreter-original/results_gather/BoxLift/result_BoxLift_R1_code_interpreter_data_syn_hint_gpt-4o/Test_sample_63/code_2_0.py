from itertools import combinations

boxes = [49, 186, 267, 243, 352, 74, 160, 115, 138, 301, 250, 145, 294, 232, 144, 293, 287, 358, 267, 266]
lifters = [60, 92, 94, 68, 72, 79]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 12:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        if box is None:
            continue
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in lifter_indices):
                    total_capacity = sum(lifters[j] for j in lifter_indices)
                    if total_capacity >= box:
                        step.append((box, list(lifter_indices)))
                        for j in lifter_indices:
                            used_lifters[j] = True
                        boxes[i] = None  # Mark this box as lifted
                        break
            if boxes[i] is None:
                break
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")