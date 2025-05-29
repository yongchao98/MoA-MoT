# List of box weights
boxes = [166, 144, 53, 213, 51, 156, 197, 311, 177, 358, 172, 134, 179, 145, 91, 188, 352, 294, 292, 88, 97, 394, 123, 294]

# List of lifters' capacities
lifters = [88, 185, 145, 195, 147, 145]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

# Function to find the best combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    from itertools import combinations
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return []

# Initialize steps
steps = []

# While there are boxes to lift
while boxes:
    step = []
    available_lifters = lifters[:]
    boxes_to_remove = []
    
    for box in boxes:
        lifters_for_box = find_lifters_for_box(box, available_lifters)
        if lifters_for_box:
            step.append((box, lifters_for_box))
            for lifter in lifters_for_box:
                available_lifters.remove(lifter)
            boxes_to_remove.append(box)
    
    # Remove lifted boxes
    for box in boxes_to_remove:
        boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")