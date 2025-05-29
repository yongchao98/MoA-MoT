from itertools import combinations

# Box weights and lifter capacities
boxes = [76, 97, 14, 30, 36, 83, 45, 46]
lifters = [69, 59, 68, 50]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Function to find combinations of lifters that can lift a given box
def find_lifter_combinations(box_weight, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return None

# Assign boxes to lifters in steps
steps = []
remaining_boxes = boxes[:]
while remaining_boxes:
    step = []
    used_lifters = set()
    for box in remaining_boxes[:]:
        lifter_indices = find_lifter_combinations(box, [l for i, l in enumerate(lifters) if i not in used_lifters])
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
            remaining_boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")