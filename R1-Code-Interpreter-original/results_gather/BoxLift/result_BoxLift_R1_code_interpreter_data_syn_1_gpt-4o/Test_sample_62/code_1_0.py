boxes = [91, 207, 152, 47, 209, 54, 251, 176, 194, 221, 152, 141, 128, 159, 57, 184]
lifters = [149, 131, 113, 109, 124]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes[:]):
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes.remove(box)
                break
        else:
            # If no single lifter can lift the box, try to combine lifters
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step.append((box, lifter_indices))
                        boxes.remove(box)
                        break
            else:
                # If we can't lift the box even with all lifters, we need to stop
                break
    steps.append(step)
    # If we have reached 6 steps, we need to ensure all boxes are lifted
    if len(steps) == 6 and boxes:
        remaining_boxes = boxes[:]
        step = []
        used_lifters = [False] * len(lifters)
        for box in remaining_boxes:
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step.append((box, lifter_indices))
                        boxes.remove(box)
                        break
        steps[-1].extend(step)
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")