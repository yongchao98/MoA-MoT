# Box weights and lifter capacities
boxes = [128, 144, 42, 49, 132, 118, 32, 48, 21, 94, 85, 51]
lifters = [85, 90, 40, 102]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter first
    for i, lifter in enumerate(available_lifters):
        if lifter >= box_weight:
            return [i]
    
    # If no single lifter can lift the box, try combinations
    for i in range(len(available_lifters)):
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box_weight:
                return [i, j]
    
    # If no combination can lift the box, return an empty list
    return []

# Initialize steps
steps = []

# Lift boxes in steps
while boxes:
    step = []
    available_lifters = lifters[:]
    lifter_indices = list(range(len(lifters)))
    
    for box in boxes[:]:
        lifter_combination = find_lifters_for_box(box, available_lifters)
        if lifter_combination:
            step.append((box, [lifter_indices[i] for i in lifter_combination]))
            # Remove used lifters from available lifters
            for i in sorted(lifter_combination, reverse=True):
                del available_lifters[i]
                del lifter_indices[i]
            # Remove the box from the list
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")