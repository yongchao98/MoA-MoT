boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
lifters = [78, 98, 51, 75, 43]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(5):
    current_step = []
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
                current_step.append((box, [j]))
                used_lifters[j] = True
                used_boxes[i] = True
                break
        else:
            # Try to find a combination of lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                total_capacity = lifters[j]
                lifter_indices = [j]
                used_lifters[j] = True
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if total_capacity + lifters[k] >= box:
                        lifter_indices.append(k)
                        used_lifters[k] = True
                        current_step.append((box, lifter_indices))
                        used_boxes[i] = True
                        break
                    total_capacity += lifters[k]
                    lifter_indices.append(k)
                    used_lifters[k] = True
                if used_boxes[i]:
                    break
    
    steps.append(current_step)
    if all(used_boxes):
        break

print(steps)