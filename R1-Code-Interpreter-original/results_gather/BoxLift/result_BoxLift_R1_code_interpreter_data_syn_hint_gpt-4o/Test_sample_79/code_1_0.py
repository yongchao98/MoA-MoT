from itertools import combinations

# Box weights and lifter capacities
boxes = [368, 220, 131, 287, 153, 175, 153, 217, 262, 115, 365, 394, 291, 222, 142, 126, 189, 82, 159, 367, 356, 168, 188, 355]
lifters = [200, 91, 128, 198, 196, 119, 188]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for indices in combinations(range(len(available_lifters)), r):
            if sum(available_lifters[i] for i in indices) >= box_weight:
                return indices
    return None

# Assign boxes to lifters in steps
steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            for i in lifter_indices:
                available_lifters[i] = 0  # Mark these lifters as used
            boxes.remove(box)
    steps.append(step)
    if len(steps) > 6:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")