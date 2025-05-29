from itertools import combinations

# Box weights and lifter capacities
boxes = [87, 123, 302, 308, 370, 86, 182, 79, 135, 170, 198, 243, 45, 212, 183, 174, 335, 134, 58, 131, 124, 360, 192, 245]
lifters = [167, 156, 74, 177, 77, 60]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find a combination of lifters that can lift a given box
def find_lifter_combination(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return None

# Assign boxes to lifters in steps
steps = []
while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        lifter_indices = find_lifter_combination(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")