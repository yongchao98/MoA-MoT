boxes = [18, 93, 23, 70, 72, 22, 96, 96]
lifters = [60, 46, 61]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    for box in boxes[:]:
        # Try to find a single lifter for the box
        for i, lifter in enumerate(remaining_lifters):
            if lifter >= box:
                step.append((box, [i]))
                remaining_lifters.pop(i)
                boxes.remove(box)
                break
        else:
            # Try to combine lifters to lift the box
            for i in range(len(remaining_lifters)):
                for j in range(i + 1, len(remaining_lifters)):
                    if remaining_lifters[i] + remaining_lifters[j] >= box:
                        step.append((box, [i, j]))
                        remaining_lifters.pop(j)
                        remaining_lifters.pop(i)
                        boxes.remove(box)
                        break
                else:
                    continue
                break
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")