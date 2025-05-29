boxes = [12, 90, 67, 10, 10, 94, 60, 46]
lifters = [56, 77, 72]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(4):
    current_step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                current_step.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # If no single lifter can lift the box, try to combine lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        current_step.append((box, [j, k]))
                        used_boxes[i] = True
                        used_lifters[j] = True
                        used_lifters[k] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(current_step)
    
    # Check if all boxes are lifted
    if all(used_boxes):
        break

# Print the steps
for step_num, step in enumerate(steps, start=1):
    print(f"Step {step_num}: {step}")