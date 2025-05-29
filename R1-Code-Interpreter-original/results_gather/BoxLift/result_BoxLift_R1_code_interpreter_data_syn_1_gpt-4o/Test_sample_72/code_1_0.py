from itertools import combinations

boxes = [280, 109, 195, 129, 54, 216, 157, 193, 44, 233, 171, 218, 204, 120, 131, 151]
lifters = [49, 131, 66, 124, 91]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    used_lifters_indices = set()
    for box in boxes[:]:
        # Try to find a combination of lifters that can lift the box
        found = False
        for r in range(1, len(available_lifters) + 1):
            for combo in combinations(enumerate(available_lifters), r):
                indices, lifter_combo = zip(*combo)
                if sum(lifter_combo) >= box and not used_lifters_indices.intersection(indices):
                    step.append((box, list(indices)))
                    used_lifters_indices.update(indices)
                    boxes.remove(box)
                    found = True
                    break
            if found:
                break
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")