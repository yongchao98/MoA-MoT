from itertools import combinations

boxes = [289, 375, 107, 145, 257, 48, 141, 83, 136, 368, 59, 133, 186, 266, 353, 73, 66, 210, 247, 79, 342, 318, 337, 162]
lifters = [75, 189, 162, 181, 118, 194, 137]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in lifter_indices):
                    if sum(lifters[i] for i in lifter_indices) >= box:
                        step.append((box, list(lifter_indices)))
                        for i in lifter_indices:
                            used_lifters[i] = True
                        boxes.remove(box)
                        break
            else:
                continue
            break
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")