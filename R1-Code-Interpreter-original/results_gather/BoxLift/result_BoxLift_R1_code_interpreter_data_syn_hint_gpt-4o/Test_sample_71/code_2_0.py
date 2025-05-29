from itertools import combinations

boxes = [156, 123, 139, 191, 91, 193, 100, 47, 36, 169, 90, 117]
lifters = [52, 55, 63, 44]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_indices = set()
    for box in boxes[:]:
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, capacities = zip(*combo)
                if sum(capacities) >= box and not used_indices.intersection(indices):
                    step.append((box, list(indices)))
                    boxes.remove(box)
                    used_indices.update(indices)
                    break
            if box not in boxes:
                break
    steps.append(step)

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")