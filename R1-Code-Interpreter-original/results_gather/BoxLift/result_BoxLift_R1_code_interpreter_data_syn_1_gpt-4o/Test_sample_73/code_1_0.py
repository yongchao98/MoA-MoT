from itertools import combinations

boxes = [75, 210, 297, 242, 276, 134, 383, 358, 146, 394, 59, 350, 94, 152, 250, 136, 116, 342, 283, 87]
lifters = [101, 62, 125, 126, 138, 54]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []

def find_lifter_combination(box_weight, lifters):
    # Try all combinations of lifters to find a valid one
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box_weight:
                return list(indices)
    return None

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    for box in boxes:
        lifter_indices = find_lifter_combination(box, [l for i, l in enumerate(lifters) if not used_lifters[i]])
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                used_lifters[index] = True
        else:
            remaining_boxes.append(box)
    steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")