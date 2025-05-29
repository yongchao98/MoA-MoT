from itertools import combinations

# Box weights and lifter capacities
boxes = [112, 109, 151, 97, 375, 276, 159, 378, 284, 275, 120, 60, 169, 398, 125, 268, 110, 92, 166, 349]
lifters = [102, 116, 199, 175, 148, 161, 161]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return combo
    return None

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
    print(f"Step {i+1}: {step}")