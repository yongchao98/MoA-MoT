# Box weights and lifter capacities
boxes = [14, 51, 14, 69, 77, 20, 52, 96]
lifters = [69, 66, 74]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Initialize steps
steps = []

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter first
    for i, capacity in enumerate(available_lifters):
        if capacity >= box_weight:
            return [i]
    
    # If no single lifter can lift the box, try combinations
    for i in range(len(available_lifters)):
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box_weight:
                return [i, j]
    
    return []

# Lift boxes in steps
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for index in lifter_indices:
                available_lifters[index] = 0  # Mark lifter as used in this step
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")