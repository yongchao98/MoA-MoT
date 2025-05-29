from itertools import combinations

boxes = [49, 239, 263, 284, 45, 159, 306, 398, 279, 137, 218, 387, 161, 287, 171, 345, 373, 271, 60, 396]
lifters = [151, 91, 70, 108, 182, 200]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_lifters_indices = set()
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, lifter_combo = zip(*combo)
                if sum(lifter_combo) >= box and not used_lifters_indices.intersection(indices):
                    step.append((box, list(indices)))
                    used_lifters_indices.update(indices)
                    boxes.remove(box)
                    break
            else:
                continue
            break
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")