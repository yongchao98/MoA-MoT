boxes = [41, 31, 60, 20, 11, 26, 52, 98]
lifters = [49, 55, 63]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(4):  # We have a maximum of 4 steps
    step_assignments = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                # Single lifter can lift the box
                step_assignments.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to find a combination of lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        step_assignments.append((box, [j, k]))
                        used_boxes[i] = True
                        used_lifters[j] = True
                        used_lifters[k] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(step_assignments)
    
    # Check if all boxes are lifted
    if all(used_boxes):
        break

print(steps)