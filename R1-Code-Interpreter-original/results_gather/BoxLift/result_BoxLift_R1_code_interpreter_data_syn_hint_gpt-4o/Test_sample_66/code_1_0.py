# Box weights and lifter capacities
boxes = [181, 52, 82, 130, 92, 195, 101, 106, 189, 42, 196, 46]
lifters = [45, 30, 52, 104, 30]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def find_min_steps(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for i, box in enumerate(boxes[:]):
            # Try to lift the box with available lifters
            lifter_indices = []
            remaining_weight = box
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and remaining_weight > 0:
                    if lifter <= remaining_weight:
                        lifter_indices.append(j)
                        remaining_weight -= lifter
                        used_lifters[j] = True
            if remaining_weight <= 0:
                step.append((box, lifter_indices))
                boxes.remove(box)
        steps.append(step)
    return steps

# Get the steps
steps = find_min_steps(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")