boxes = [12, 90, 67, 10, 10, 94, 60, 46]
lifters = [56, 77, 72]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(4):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lift.append((box, [j]))
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
                        step_lift.append((box, [j, k]))
                        used_boxes[i] = True
                        used_lifters[j] = True
                        used_lifters[k] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(step_lift)
    
    if all(used_boxes):
        break

print(steps)