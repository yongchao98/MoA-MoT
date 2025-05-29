boxes = [280, 32, 168, 120, 112, 285, 46, 33, 148, 74, 295, 279, 207, 138, 153, 63]
lifters = [53, 92, 104, 86, 106, 68]

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
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and remaining_weight > 0:
                if lifter <= remaining_weight:
                    lifter_indices.append(j)
                    remaining_weight -= lifter
                    used_lifters[j] = True
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes[i] = None
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")