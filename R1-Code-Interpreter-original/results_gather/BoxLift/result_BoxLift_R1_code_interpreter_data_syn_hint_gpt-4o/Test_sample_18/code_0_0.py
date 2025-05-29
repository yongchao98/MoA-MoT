boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
lifters = [122, 116, 51, 165, 95, 150, 179]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if not boxes:
            break
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
                for k, lifter2 in enumerate(lifters):
                    if not used_lifters[k]:
                        total_capacity += lifter2
                        lifter_indices.append(k)
                        if total_capacity >= box:
                            for idx in lifter_indices:
                                used_lifters[idx] = True
                            step.append((box, lifter_indices))
                            boxes[i] = None
                            break
                if boxes[i] is None:
                    break
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")