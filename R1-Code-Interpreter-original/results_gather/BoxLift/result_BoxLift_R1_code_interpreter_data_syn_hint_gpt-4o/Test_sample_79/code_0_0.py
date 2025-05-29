# Box weights and lifter capacities
boxes = [368, 220, 131, 287, 153, 175, 153, 217, 262, 115, 365, 394, 291, 222, 142, 126, 189, 82, 159, 367, 356, 168, 188, 355]
lifters = [200, 91, 128, 198, 196, 119, 188]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    used_indices = []
    total_capacity = 0
    for i, capacity in enumerate(available_lifters):
        if total_capacity >= box_weight:
            break
        if capacity > 0:
            used_indices.append(i)
            total_capacity += capacity
            available_lifters[i] = 0  # Mark this lifter as used
    if total_capacity >= box_weight:
        return used_indices
    return None

# Assign boxes to lifters in steps
steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
    steps.append(step)
    if len(steps) > 6:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")