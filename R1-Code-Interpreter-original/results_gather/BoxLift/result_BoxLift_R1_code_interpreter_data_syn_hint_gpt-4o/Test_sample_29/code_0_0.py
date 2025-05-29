boxes = [100, 225, 68, 49, 127, 244, 299, 344, 330, 310, 299, 189, 124, 216, 350, 367, 85, 231, 165, 144, 323, 342, 368, 179]
lifters = [64, 93, 158, 137, 107, 95, 173]

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