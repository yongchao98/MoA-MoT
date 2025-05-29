# Box weights and lifter capacities
boxes = [64, 70, 386, 351, 113, 77, 314, 333, 266, 399, 193, 44, 181, 200, 238, 175, 370, 118, 337, 134]
lifters = [140, 115, 159, 147, 129, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, lifters):
    from itertools import combinations
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return []

# Assign boxes to lifters
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