# Box weights and lifter capacities
boxes = [18, 27, 65, 31, 60, 73, 30, 67]
lifters = [66, 72, 43]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter who can lift the box
    for i, capacity in enumerate(available_lifters):
        if capacity >= box_weight:
            return [i]
    
    # If no single lifter can lift the box, try combinations
    for i in range(len(available_lifters)):
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box_weight:
                return [i, j]
    
    return None

# Initialize steps
steps = []

# Try to lift all boxes in 4 steps or less
for step in range(4):
    step_lift = []
    available_lifters = lifters[:]
    
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step_lift.append((box, lifter_indices))
            boxes.remove(box)
            for index in lifter_indices:
                available_lifters[index] = 0  # Mark lifter as used for this step
    
    steps.append(step_lift)
    
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")