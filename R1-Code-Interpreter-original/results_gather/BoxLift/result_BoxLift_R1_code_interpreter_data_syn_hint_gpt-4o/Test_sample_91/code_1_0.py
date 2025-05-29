from itertools import combinations

# Box weights and lifter capacities
boxes = [159, 180, 64, 184, 110, 43, 128, 47, 47, 138, 62, 97]
lifters = [81, 98, 62, 46]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Function to find the combination of lifters that can lift a given box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box_weight:
                return list(indices)
    return None

# Plan the steps
steps = []
while boxes and len(steps) < 7:
    step = []
    available_lifters = lifters[:]
    used_indices = set()
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for index in lifter_indices:
                used_indices.add(index)
    # Mark used lifters as unavailable for this step
    for index in used_indices:
        available_lifters[index] = 0
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")