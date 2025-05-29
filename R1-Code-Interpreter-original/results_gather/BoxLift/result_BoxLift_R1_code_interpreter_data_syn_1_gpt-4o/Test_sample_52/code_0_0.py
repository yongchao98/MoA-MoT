# Box weights and lifter capacities
boxes = [75, 298, 70, 200, 278, 35, 278, 56, 270, 212, 161, 98, 244, 292, 177, 145]
lifters = [105, 102, 125, 93, 123]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    from itertools import combinations
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return []

# Assign boxes to lifters in steps
steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_combo = find_lifters_for_box(box, available_lifters)
        if lifter_combo:
            step.append((box, [lifters.index(l) for l in lifter_combo]))
            boxes.remove(box)
            for l in lifter_combo:
                available_lifters.remove(l)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")