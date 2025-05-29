from itertools import combinations

# Box weights and lifter capacities
boxes = [44, 178, 94, 176, 159, 199, 197, 133, 173, 186, 196, 80]
lifters = [61, 45, 66, 38, 52]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a combination of lifters that can lift a given box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return None

steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_combo = find_lifters_for_box(box, available_lifters)
        if lifter_combo:
            step.append((box, lifter_combo))
            boxes.remove(box)
            for lifter in lifter_combo:
                available_lifters.remove(lifter)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")