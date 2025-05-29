boxes = [18, 27, 65, 31, 60, 73, 30, 67]
lifters = [66, 72, 43]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(4):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lifting.append((box, [j]))
                used_lifters[j] = True
                used_boxes[i] = True
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
                        step_lifting.append((box, [j, k]))
                        used_lifters[j] = True
                        used_lifters[k] = True
                        used_boxes[i] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(step_lifting)
    if all(used_boxes):
        break

print(steps)