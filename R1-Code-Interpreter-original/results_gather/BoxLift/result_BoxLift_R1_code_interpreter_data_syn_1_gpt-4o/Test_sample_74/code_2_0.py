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
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo) and sum(lifters[j] for j in combo) >= box:
                    step.append((box, list(combo)))
                    for j in combo:
                        used_lifters[j] = True
                    boxes[i] = None  # Mark this box as lifted
                    break
            if boxes[i] is None:
                break
    boxes = [b for b in boxes if b is not None]  # Remove lifted boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")