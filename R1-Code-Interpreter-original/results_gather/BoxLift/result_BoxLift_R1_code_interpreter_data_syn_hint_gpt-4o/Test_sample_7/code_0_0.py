boxes = [202, 297, 200, 241, 156, 177, 70, 78, 133, 176, 75, 70, 167, 240, 299, 91]
lifters = [56, 160, 118, 138, 160, 150]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(5):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if used_lifters[j]:
                continue
            if lifter >= box:
                # Single lifter can lift the box
                step_lifting.append((box, [j]))
                used_lifters[j] = True
                used_boxes[i] = True
                break
        else:
            # Try to find a combination of lifters
            lifter_indices = []
            total_capacity = 0
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    lifter_indices.append(j)
                    total_capacity += lifter
                    if total_capacity >= box:
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        step_lifting.append((box, lifter_indices))
                        used_boxes[i] = True
                        break
    
    steps.append(step_lifting)
    if all(used_boxes):
        break

print(steps)