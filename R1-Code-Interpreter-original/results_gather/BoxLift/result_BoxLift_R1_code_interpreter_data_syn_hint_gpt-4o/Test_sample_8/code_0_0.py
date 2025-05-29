# Box weights and lifter capacities
boxes = [362, 202, 192, 149, 357, 335, 209, 193, 334, 375, 374, 167, 276, 70, 282, 201, 242, 158, 236, 212]
lifters = [174, 130, 182, 186, 80, 149, 87]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the optimal assignment of lifters to boxes
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for i, box in enumerate(boxes[:]):
            # Try to find a single lifter or a combination of lifters to lift the box
            lifter_indices = []
            remaining_weight = box
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter <= remaining_weight:
                    lifter_indices.append(j)
                    remaining_weight -= lifter
                    used_lifters[j] = True
                if remaining_weight <= 0:
                    break
            if remaining_weight <= 0:
                step.append((box, lifter_indices))
                boxes.remove(box)
        steps.append(step)
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")