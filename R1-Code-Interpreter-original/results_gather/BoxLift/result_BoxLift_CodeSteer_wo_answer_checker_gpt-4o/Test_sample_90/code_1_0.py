boxes = [71, 55, 11, 10, 71, 91, 14, 69]
lifters = [46, 54, 47]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
max_steps = 5

while boxes and len(steps) < max_steps:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = boxes[:]
    
    for i, box in enumerate(remaining_boxes):
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
if boxes:
    output = "<<<No valid solution found within the constraints>>>"
else:
    output = "<<<"
    for i, step in enumerate(steps):
        output += f"Step {i + 1}: {step}\n"
    output += ">>>"

print(output)