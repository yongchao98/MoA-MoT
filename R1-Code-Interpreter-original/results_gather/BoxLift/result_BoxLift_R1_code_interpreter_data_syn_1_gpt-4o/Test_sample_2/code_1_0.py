boxes = [66, 299, 90, 133, 52, 245, 57, 62, 71, 217, 117, 193, 204, 84, 224, 51]
lifters = [123, 155, 65, 92, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    lifted_indices = []

    for i, box in enumerate(boxes):
        if i in lifted_indices:
            continue
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                lifted_indices.append(i)
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
                    if total_capacity + lifters[k] >= box:
                        lifter_indices.append(k)
                        step.append((box, lifter_indices))
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        lifted_indices.append(i)
                        break
                    total_capacity += lifters[k]
                    lifter_indices.append(k)
                if i in lifted_indices:
                    break

    # Remove lifted boxes by index
    boxes = [box for i, box in enumerate(boxes) if i not in lifted_indices]
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")