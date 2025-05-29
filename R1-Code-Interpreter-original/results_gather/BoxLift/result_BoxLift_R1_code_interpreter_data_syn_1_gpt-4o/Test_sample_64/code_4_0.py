from itertools import combinations

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
        if box is None:
            continue
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo) and sum(lifters[j] for j in combo) >= box:
                    for j in combo:
                        used_lifters[j] = True
                    step.append((box, list(combo)))
                    boxes[i] = None  # Mark the box as lifted
                    break
            if boxes[i] is None:
                break
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")