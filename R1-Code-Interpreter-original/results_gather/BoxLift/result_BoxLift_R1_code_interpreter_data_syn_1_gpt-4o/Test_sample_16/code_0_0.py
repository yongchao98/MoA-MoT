# Box weights and lifter capacities
boxes = [166, 194, 27, 144, 198, 172, 81, 117, 154, 183, 27, 184]
lifters = [97, 90, 80, 111, 76]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for i, box in enumerate(boxes[:]):
            # Try to find a single lifter or combination of lifters to lift the box
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    boxes.remove(box)
                    break
            else:
                # Try combinations of lifters
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
        if len(steps) > 5:
            return "Not possible to lift all boxes in 5 steps"
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
print(steps)