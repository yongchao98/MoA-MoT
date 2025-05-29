from itertools import combinations

# Box weights and lifter capacities
boxes = [254, 205, 272, 360, 320, 328, 267, 223, 337, 238, 134, 235, 388, 160, 244, 378, 349, 262, 104, 396]
lifters = [132, 86, 81, 199, 93, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of lifters needed to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return []

# Assign boxes to lifters
steps = []
while boxes:
    step = []
    used_lifters = []
    for box in boxes[:]:
        lifters_for_box = find_lifters_for_box(box, [l for l in lifters if l not in used_lifters])
        if lifters_for_box:
            step.append((box, [lifters.index(l) for l in lifters_for_box]))
            used_lifters.extend(lifters_for_box)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")