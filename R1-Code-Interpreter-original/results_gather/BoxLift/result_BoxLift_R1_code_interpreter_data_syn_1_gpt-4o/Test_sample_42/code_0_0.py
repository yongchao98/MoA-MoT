boxes = [86, 183, 78, 85, 112, 170, 130, 121, 208, 226, 39, 176, 68, 256, 56, 34]
lifters = [116, 81, 149, 77, 138]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

while remaining_boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to find a single lifter or a combination of lifters to lift the box
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and lifter >= box:
                step.append((box, [i]))
                used_lifters[i] = True
                remaining_boxes.remove(box)
                break
        else:
            # Try to find a combination of lifters
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
                        remaining_boxes.remove(box)
                        break
                if box not in remaining_boxes:
                    break
    
    steps.append(step)

print(steps)