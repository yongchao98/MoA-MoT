boxes = [11, 87, 70, 21, 90, 54, 50, 27, 13, 76]
lifters = [72, 52, 75]

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
            # Try to combine lifters if a single one can't lift the box
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

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")