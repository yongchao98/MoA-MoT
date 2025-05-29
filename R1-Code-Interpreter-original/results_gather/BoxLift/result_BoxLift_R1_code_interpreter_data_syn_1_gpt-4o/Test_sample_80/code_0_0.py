boxes = [14, 13, 59, 47, 73, 64, 58, 75]
lifters = [74, 57, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for _ in range(4):  # Maximum 4 steps
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to combine lifters to lift the box
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        step.append((box, [j, k]))
                        used_boxes[i] = True
                        used_lifters[j] = True
                        used_lifters[k] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(step)
    if all(used_boxes):
        break

print(steps)