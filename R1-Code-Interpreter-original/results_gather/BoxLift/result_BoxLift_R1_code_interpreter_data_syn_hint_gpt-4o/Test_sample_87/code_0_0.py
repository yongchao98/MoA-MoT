boxes = [169, 198, 165, 154, 71, 159, 55, 205, 299, 170, 122, 160, 43, 259, 246, 172]
lifters = [80, 51, 45, 141, 83, 152]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(7):
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
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to find a combination of lifters
            lifter_indices = []
            total_capacity = 0
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    if total_capacity >= box:
                        step_lifting.append((box, lifter_indices))
                        used_boxes[i] = True
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        break
    
    steps.append(step_lifting)
    if all(used_boxes):
        break

# Print the steps
for step_num, step in enumerate(steps, start=1):
    print(f"Step {step_num}: {step}")