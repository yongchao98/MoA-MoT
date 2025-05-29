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
        for i, box in enumerate(boxes):
            # Try to find a single lifter or combination of lifters to lift the box
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    break
            else:
                # Try combinations if a single lifter cannot lift the box
                for j in range(len(lifters)):
                    if used_lifters[j]:
                        continue
                    total_capacity = lifters[j]
                    lifter_indices = [j]
                    for k in range(j + 1, len(lifters)):
                        if used_lifters[k]:
                            continue
                        if total_capacity + lifters[k] >= box:
                            lifter_indices.append(k)
                            step.append((box, lifter_indices))
                            for idx in lifter_indices:
                                used_lifters[idx] = True
                            break
                        total_capacity += lifters[k]
                        lifter_indices.append(k)
                    if len(step) > len(boxes) - i:
                        break
        # Remove lifted boxes
        boxes = [box for box, _ in step if box not in [b for b, _ in step]]
        steps.append(step)
        if len(steps) > 4:
            break
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")