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
        # Try to combine lifters to lift the box
        total_capacity = 0
        lifter_indices = []
        for j, lifter in enumerate(lifters):
            if not used_lifters[j]:
                total_capacity += lifter
                lifter_indices.append(j)
                if total_capacity >= box:
                    for idx in lifter_indices:
                        used_lifters[idx] = True
                    step.append((box, lifter_indices))
                    boxes[i] = None
                    break
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")