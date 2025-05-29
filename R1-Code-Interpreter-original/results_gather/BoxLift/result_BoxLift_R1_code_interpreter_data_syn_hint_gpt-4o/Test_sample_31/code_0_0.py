boxes = [44, 70, 139, 184, 194, 67, 188, 33, 149, 194, 132, 117]
lifters = [110, 109, 30, 34]

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
            # Try combinations of lifters
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

print(steps)