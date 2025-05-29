boxes = [16, 56, 15, 82, 85, 12, 12, 32]
lifters = [63, 48, 72, 68]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(3):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lifting.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # If no single lifter can lift the box, try to combine lifters
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step_lifting.append((box, lifter_indices))
                        used_boxes[i] = True
                        break
    
    steps.append(step_lifting)
    if all(used_boxes):
        break

print(steps)