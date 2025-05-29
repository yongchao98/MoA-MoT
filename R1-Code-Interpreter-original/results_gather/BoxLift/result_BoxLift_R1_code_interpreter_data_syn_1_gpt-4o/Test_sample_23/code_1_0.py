from itertools import combinations

boxes = [254, 205, 272, 360, 320, 328, 267, 223, 337, 238, 134, 235, 388, 160, 244, 378, 349, 262, 104, 396]
lifters = [132, 86, 81, 199, 93, 92]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        # Find a combination of lifters to lift the box
        lifter_indices = []
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[i] for i in combo):
                    if sum(lifters[i] for i in combo) >= box:
                        lifter_indices = list(combo)
                        break
            if lifter_indices:
                break
        if lifter_indices:
            step.append((box, lifter_indices))
            for i in lifter_indices:
                used_lifters[i] = True
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")