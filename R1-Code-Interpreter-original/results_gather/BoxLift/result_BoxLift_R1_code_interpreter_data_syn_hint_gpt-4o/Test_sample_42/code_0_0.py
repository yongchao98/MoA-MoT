boxes = [86, 183, 78, 85, 112, 170, 130, 121, 208, 226, 39, 176, 68, 256, 56, 34]
lifters = [116, 81, 149, 77, 138]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(6):
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
            lifter_indices = []
            total_capacity = 0
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    lifter_indices.append(j)
                    total_capacity += lifter
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step_lifting.append((box, lifter_indices))
                        used_boxes[i] = True
                        break
            else:
                # If we can't lift this box in this step, reset the lifters used
                for idx in lifter_indices:
                    used_lifters[idx] = False
    
    if not step_lifting:
        break
    
    steps.append(step_lifting)

# Check if all boxes are lifted
if all(used_boxes):
    for idx, step in enumerate(steps):
        print(f"Step {idx + 1}: {step}")
else:
    print("Not all boxes can be lifted in 6 steps.")