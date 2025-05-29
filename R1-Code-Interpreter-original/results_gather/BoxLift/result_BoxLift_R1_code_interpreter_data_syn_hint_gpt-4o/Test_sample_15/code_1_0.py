from itertools import combinations

# Box weights and lifter capacities
boxes = [50, 70, 83, 17, 21, 27, 53, 31]
lifters = [42, 53, 74, 78]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find lifter combinations that can lift a given box
def find_lifter_combination(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box_weight:
                return list(indices)
    return None

# Assign boxes to lifters in steps
steps = []
remaining_boxes = boxes[:]
while remaining_boxes and len(steps) < 3:
    step = []
    used_lifters = set()
    for box in remaining_boxes[:]:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        lifter_indices = find_lifter_combination(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
            remaining_boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")