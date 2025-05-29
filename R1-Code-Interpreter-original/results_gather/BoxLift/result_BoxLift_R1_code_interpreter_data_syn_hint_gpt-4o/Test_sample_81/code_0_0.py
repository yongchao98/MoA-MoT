from itertools import combinations

# Box weights and lifter capacities
boxes = [174, 117, 74, 120, 98, 65, 112, 61, 46, 87, 176, 110]
lifters = [119, 104, 67, 98]

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
    remaining_boxes = []
    
    for box in boxes:
        lifter_combo = find_lifters_for_box(box, available_lifters)
        if lifter_combo:
            step.append((box, [lifters.index(l) for l in lifter_combo]))
            for l in lifter_combo:
                available_lifters.remove(l)
        else:
            remaining_boxes.append(box)
    
    steps.append(step)
    boxes = remaining_boxes
    
    if len(steps) > 5:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")