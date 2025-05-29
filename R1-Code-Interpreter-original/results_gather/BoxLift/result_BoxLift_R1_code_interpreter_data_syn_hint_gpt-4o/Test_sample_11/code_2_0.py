from itertools import combinations

# Box weights and lifter capacities
boxes = [196, 90, 142, 22, 143, 101, 43, 166, 79, 129, 77, 187]
lifters = [71, 58, 104, 105]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
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
            step.append((box, lifter_combo))
            boxes.remove(box)
            for lifter in lifter_combo:
                available_lifters.remove(lifter)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")