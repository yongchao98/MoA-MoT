from itertools import combinations

# Box weights and lifter capacities
boxes = [260, 319, 78, 313, 164, 297, 367, 392, 175, 245, 200, 45, 386, 124, 169, 97, 359, 273, 175, 238]
lifters = [57, 82, 155, 87, 174, 192, 149]

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

steps = []
while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        lifter_combo = find_lifters_for_box(box, available_lifters)
        if lifter_combo:
            step.append((box, [lifters.index(l) for l in lifter_combo]))
            used_lifters.update(lifters.index(l) for l in lifter_combo)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")