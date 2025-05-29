from itertools import combinations

# Box weights and lifter capacities
boxes = [41, 69, 80, 19, 54, 85, 64, 35]
lifters = [45, 45, 62, 74]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        remaining_boxes = []
        for box in boxes:
            # Try to find a single lifter or combination of lifters to lift the box
            lifter_indices = []
            total_capacity = 0
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    if total_capacity >= box:
                        step.append((box, lifter_indices))
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        break
            else:
                # If no combination can lift the box, add it to remaining boxes
                remaining_boxes.append(box)
        boxes = remaining_boxes
        steps.append(step)
        if len(steps) > 4:
            break
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")