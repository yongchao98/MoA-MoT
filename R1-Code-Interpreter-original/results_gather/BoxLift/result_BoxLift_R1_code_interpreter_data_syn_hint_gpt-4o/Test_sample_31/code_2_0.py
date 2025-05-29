boxes = [44, 70, 139, 184, 194, 67, 188, 33, 149, 194, 132, 117]
lifters = [110, 109, 30, 34]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

steps = []

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a single lifter or a combination of lifters to lift the box
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and lifter >= box:
                step.append((box, [i]))
                used_lifters[i] = True
                boxes.remove(box)
                break
        else:
            # Try combinations of lifters
            for i in range(len(lifters)):
                if used_lifters[i]:
                    continue
                for j in range(i + 1, len(lifters)):
                    if used_lifters[j]:
                        continue
                    if lifters[i] + lifters[j] >= box:
                        step.append((box, [i, j]))
                        used_lifters[i] = True
                        used_lifters[j] = True
                        boxes.remove(box)
                        break
                else:
                    continue
                break
    
    steps.append(step)

print(steps)