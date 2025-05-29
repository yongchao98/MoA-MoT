boxes = [112, 109, 151, 97, 375, 276, 159, 378, 284, 275, 120, 60, 169, 398, 125, 268, 110, 92, 166, 349]
lifters = [102, 116, 199, 175, 148, 161, 161]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes[:]):
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes.remove(box)
                break
        else:
            # Try to combine lifters if a single lifter can't lift the box
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
                        boxes.remove(box)
                        break
                    total_capacity += lifters[k]
                    lifter_indices.append(k)
                if box not in boxes:
                    break
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")