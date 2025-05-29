# Define the weights of the boxes and the capacities of the lifters
boxes = [21, 96, 83, 73, 100, 73, 55, 25]
lifters = [62, 79, 54]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Initialize the steps
steps = []

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter or a combination of lifters to lift the box
    for i in range(len(available_lifters)):
        if available_lifters[i] >= box_weight:
            return [i]
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box_weight:
                return [i, j]
    return []

# Lift boxes in steps
while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, [lifters[i] for i in range(len(lifters)) if i not in used_lifters])
        if lifter_indices:
            step.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")