boxes = [38, 74, 55, 11, 82, 43, 34, 20]
lifters = [76, 78, 46, 63]

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
    
    # Check if all boxes are lifted
    if all(used_boxes):
        break

print(steps)