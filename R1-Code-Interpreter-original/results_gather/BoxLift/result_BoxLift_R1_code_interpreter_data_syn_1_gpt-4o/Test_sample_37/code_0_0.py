boxes = [97, 164, 187, 145, 255, 257, 240, 55, 261, 185, 252, 257, 206, 298, 121, 202]
lifters = [131, 106, 144, 121, 109, 52]

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
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes[i] = None
                break
        else:
            # If no single lifter can lift the box, try to combine lifters
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step.append((box, lifter_indices))
                        boxes[i] = None
                        break
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")