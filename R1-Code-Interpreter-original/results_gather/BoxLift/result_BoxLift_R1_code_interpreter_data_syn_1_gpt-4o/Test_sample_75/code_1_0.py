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
while remaining_boxes and len(steps) < 3:
    step = []
    used_lifters = set()
    for box in remaining_boxes[:]:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        lifter_indices = find_lifter_combinations(box, available_lifters)
        if lifter_indices is not None:
            # Map back to original lifter indices
            original_indices = [i for i, l in enumerate(lifters) if l in available_lifters]
            mapped_indices = [original_indices[i] for i in lifter_indices]
            step.append((box, mapped_indices))
            used_lifters.update(mapped_indices)
            remaining_boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")