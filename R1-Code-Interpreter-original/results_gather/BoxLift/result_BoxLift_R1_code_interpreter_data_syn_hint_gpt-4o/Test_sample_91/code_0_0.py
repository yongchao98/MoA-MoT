# Box weights and lifter capacities
boxes = [159, 180, 64, 184, 110, 43, 128, 47, 47, 138, 62, 97]
lifters = [81, 98, 62, 46]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Function to find the combination of lifters that can lift a given box
def find_lifters_for_box(box_weight, available_lifters):
    used_lifters = []
    remaining_weight = box_weight
    for i, lifter in enumerate(available_lifters):
        if remaining_weight <= 0:
            break
        if lifter <= remaining_weight:
            used_lifters.append(i)
            remaining_weight -= lifter
    if remaining_weight > 0:
        return None  # Not possible to lift this box with available lifters
    return used_lifters

# Plan the steps
steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for index in lifter_indices:
                available_lifters[index] = 0  # Mark lifter as used
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")