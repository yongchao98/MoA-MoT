boxes = [112, 109, 151, 97, 375, 276, 159, 378, 284, 275, 120, 60, 169, 398, 125, 268, 110, 92, 166, 349]
lifters = [102, 116, 199, 175, 148, 161, 161]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                break
        else:
            # Try to find a combination of lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                total_capacity = lifters[j]
                lifter_indices = [j]
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if total_capacity >= box:
                        break
                    total_capacity += lifters[k]
                    lifter_indices.append(k)
                if total_capacity >= box:
                    for idx in lifter_indices:
                        used_lifters[idx] = True
                    step.append((box, lifter_indices))
                    break
    # Remove lifted boxes
    boxes = [box for box, _ in step if box not in [b for b, _ in step]]
    steps.append(step)
    if len(steps) > 6:
        break

print(steps)