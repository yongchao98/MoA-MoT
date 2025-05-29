# Box weights and lifter capacities
boxes = [66, 299, 90, 133, 52, 245, 57, 62, 71, 217, 117, 193, 204, 84, 224, 51]
lifters = [123, 155, 65, 92, 92]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, lifters):
    from itertools import combinations
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return []

# Assign boxes to lifters in steps
steps = []
while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, [l for i, l in enumerate(lifters) if i not in used_lifters])
        if lifter_indices:
            step.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")