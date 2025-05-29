boxes = [66, 299, 90, 133, 52, 245, 57, 62, 71, 217, 117, 193, 204, 84, 224, 51]
lifters = [123, 155, 65, 92, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a single lifter or a combination of lifters to lift the box
        for j in range(len(available_lifters)):
            if used_lifters[j]:
                continue
            if available_lifters[j] >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes.remove(box)
                break
        else:
            # Try to combine lifters
            for j in range(len(available_lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(available_lifters)):
                    if used_lifters[k]:
                        continue
                    if available_lifters[j] + available_lifters[k] >= box:
                        step.append((box, [j, k]))
                        used_lifters[j] = True
                        used_lifters[k] = True
                        boxes.remove(box)
                        break
                if used_lifters[j]:
                    break

    steps.append(step)
    if len(steps) > 6:
        break

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")