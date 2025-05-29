# Box weights and lifter capacities
boxes = [58, 45, 94, 26, 49, 153, 61, 23, 170, 143, 75, 166]
lifters = [57, 61, 104, 98, 70]

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
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    boxes.remove(box)
                    break
            else:
                # Try combinations of lifters if a single lifter cannot lift the box
                for j in range(len(lifters)):
                    if used_lifters[j]:
                        continue
                    for k in range(j + 1, len(lifters)):
                        if used_lifters[k]:
                            continue
                        if lifters[j] + lifters[k] >= box:
                            step.append((box, [j, k]))
                            used_lifters[j] = True
                            used_lifters[k] = True
                            boxes.remove(box)
                            break
                    if used_lifters[j]:
                        break
        steps.append(step)
        if len(steps) > 4:
            break
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")