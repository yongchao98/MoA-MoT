# Box weights and lifter capacities
boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
lifters = [93, 71, 82, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    from itertools import combinations
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return []

# Assign lifters to boxes in steps
steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for lifter in lifter_indices:
                available_lifters.remove(lifter)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")