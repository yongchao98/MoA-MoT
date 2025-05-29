boxes = [49, 239, 263, 284, 45, 159, 306, 398, 279, 137, 218, 387, 161, 287, 171, 345, 373, 271, 60, 396]
lifters = [151, 91, 70, 108, 182, 200]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    for box in boxes[:]:
        # Try to find a single lifter or a combination of lifters to lift the box
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