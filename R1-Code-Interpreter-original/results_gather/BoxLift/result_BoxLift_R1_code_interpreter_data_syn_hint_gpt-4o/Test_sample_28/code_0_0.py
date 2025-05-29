boxes = [71, 72, 167, 64, 339, 45, 243, 53, 52, 280, 164, 322, 362, 75, 158, 142, 150, 353, 182, 282]
lifters = [139, 55, 172, 82, 129, 59, 130]

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
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lifting.append((box, [j]))
                used_lifters[j] = True
                used_boxes[i] = True
                break
        else:
            # Try to combine lifters for the box
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    if total_capacity >= box:
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        step_lifting.append((box, lifter_indices))
                        used_boxes[i] = True
                        break
    
    if not step_lifting:
        break
    
    steps.append(step_lifting)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")