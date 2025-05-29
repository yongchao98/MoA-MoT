from itertools import combinations

boxes = [49, 239, 263, 284, 45, 159, 306, 398, 279, 137, 218, 387, 161, 287, 171, 345, 373, 271, 60, 396]
lifters = [151, 91, 70, 108, 182, 200]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

steps = []

while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        # Try to find the smallest combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if sum(lifters[i] for i in combo) >= box and not used_lifters.intersection(combo):
                    step.append((box, list(combo)))
                    used_lifters.update(combo)
                    boxes.remove(box)
                    break
            else:
                continue
            break
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")