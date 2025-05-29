from itertools import combinations

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
    # Check single lifters first
    for i, lifter in enumerate(available_lifters):
        if lifter >= box_weight:
            return [i]
    # Check combinations of two lifters
    for i, j in combinations(range(len(available_lifters)), 2):
        if available_lifters[i] + available_lifters[j] >= box_weight:
            return [i, j]
    return []

# Try to lift all boxes in 4 steps or less
for step in range(4):
    step_assignments = []
    used_lifters = set()
    remaining_boxes = []
    for box in boxes:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices:
            step_assignments.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
        else:
            remaining_boxes.append(box)
    steps.append(step_assignments)
    boxes = remaining_boxes
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")