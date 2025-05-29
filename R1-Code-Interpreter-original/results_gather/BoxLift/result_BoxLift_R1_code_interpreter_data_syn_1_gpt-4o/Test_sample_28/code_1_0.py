from itertools import combinations

boxes = [71, 72, 167, 64, 339, 45, 243, 53, 52, 280, 164, 322, 362, 75, 158, 142, 150, 353, 182, 282]
lifters = [139, 55, 172, 82, 129, 59, 130]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[j] for j in combo):
                    if sum(lifters[j] for j in combo) >= box:
                        lifter_indices = list(combo)
                        break
            if lifter_indices:
                break
        if lifter_indices:
            for j in lifter_indices:
                used_lifters[j] = True
            step.append((box, lifter_indices))
            boxes[i] = None  # Mark this box as lifted
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")