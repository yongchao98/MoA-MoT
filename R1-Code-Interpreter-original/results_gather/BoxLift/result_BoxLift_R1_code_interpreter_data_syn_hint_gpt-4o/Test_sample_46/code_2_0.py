# Box weights and lifter capacities
boxes = [45, 32, 84, 101, 161, 118, 30, 152, 159, 20, 110, 103]
lifters = [39, 32, 112, 102, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters to lift a box
def find_lifter_combination(box_weight, lifters):
    from itertools import combinations
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return None

# Assign boxes to lifters in steps
steps = []
while boxes and len(steps) < 5:
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