# Box weights and lifter capacities
boxes = [227, 106, 388, 173, 63, 178, 261, 205, 323, 124, 386, 81, 71, 127, 134, 212, 150, 114, 41, 277, 123, 152, 47, 47]
lifters = [166, 106, 122, 181, 102, 119, 147]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    from itertools import combinations
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
    i = 0
    while i < len(boxes):
        box_weight = boxes[i]
        lifter_combo = find_lifters_for_box(box_weight, available_lifters)
        if lifter_combo:
            step.append((box_weight, [lifters.index(l) for l in lifter_combo]))
            for l in lifter_combo:
                available_lifters.remove(l)
            boxes.pop(i)
        else:
            i += 1
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")