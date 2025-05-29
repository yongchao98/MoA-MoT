from itertools import combinations

boxes = [169, 106, 153, 20, 138, 80, 94, 190, 142, 196, 186, 87]
lifters = [55, 85, 43, 100, 47]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

def find_combination(box_weight, lifters):
    # Try all combinations of lifters to find a valid one
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return None

for _ in range(7):  # Maximum 7 steps
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    for box_weight in boxes:
        lifter_indices = find_combination(box_weight, [c for i, c in enumerate(lifters) if not used_lifters[i]])
        if lifter_indices is not None:
            step.append((box_weight, lifter_indices))
            for idx in lifter_indices:
                used_lifters[idx] = True
        else:
            remaining_boxes.append(box_weight)
    steps.append(step)
    boxes = remaining_boxes
    if not boxes:
        break

print(steps)