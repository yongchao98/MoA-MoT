from itertools import combinations

boxes = [364, 138, 172, 200, 44, 194, 195, 282, 189, 400, 233, 307, 71, 300, 306, 346, 61, 136, 168, 392, 97, 367, 118, 354]
lifters = [54, 124, 186, 103, 122, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 10:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        if box is None:
            continue
        lifter_indices = []
        remaining_weight = box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo) and sum(lifters[j] for j in combo) >= remaining_weight:
                    lifter_indices = list(combo)
                    for j in combo:
                        used_lifters[j] = True
                    break
            if lifter_indices:
                break
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes[i] = None
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")