from itertools import combinations

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