boxes = [88, 313, 344, 387, 392, 370, 99, 389, 340, 391, 131, 147, 334, 251, 157, 144, 290, 167, 141, 79]
lifters = [148, 193, 54, 146, 106, 53, 88]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if not boxes:
            break
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if used_lifters[j]:
                continue
            if lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes[i] = None
                break
        else:
            # Try to combine lifters
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if used_lifters[j]:
                    continue
                if total_capacity + lifter >= box:
                    lifter_indices.append(j)
                    step.append((box, lifter_indices))
                    for idx in lifter_indices:
                        used_lifters[idx] = True
                    boxes[i] = None
                    break
                else:
                    total_capacity += lifter
                    lifter_indices.append(j)
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")