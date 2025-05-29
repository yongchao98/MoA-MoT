# Define the weights of the boxes and the capacities of the lifters
boxes = [14, 51, 14, 69, 77, 20, 52, 96]
lifters = [69, 66, 74]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Initialize steps
steps = []

# Function to find a combination of lifters that can lift a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter or a combination of lifters that can lift the box
    for i in range(len(available_lifters)):
        if available_lifters[i] >= box_weight:
            return [i], available_lifters[i]
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box_weight:
                return [i, j], available_lifters[i] + available_lifters[j]
    return [], 0

# Try to lift all boxes in 4 steps or less
for step in range(4):
    step_assignments = []
    used_lifters = set()
    for box in boxes:
        lifter_indices, total_capacity = find_lifters_for_box(box, [l for i, l in enumerate(lifters) if i not in used_lifters])
        if lifter_indices:
            step_assignments.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
    steps.append(step_assignments)
    # Remove lifted boxes
    boxes = [box for box, _ in step_assignments if box not in boxes]
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")