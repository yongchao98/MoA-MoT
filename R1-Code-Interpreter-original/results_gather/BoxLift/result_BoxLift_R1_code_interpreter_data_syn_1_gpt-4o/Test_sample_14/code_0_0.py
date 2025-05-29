boxes = [364, 138, 172, 200, 44, 194, 195, 282, 189, 400, 233, 307, 71, 300, 306, 346, 61, 136, 168, 392, 97, 367, 118, 354]
lifters = [54, 124, 186, 103, 122, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                remaining_weight -= lifter
                used_lifters[j] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
    
    # Remove lifted boxes
    for box, _ in step:
        boxes.remove(box)
    
    steps.append(step)
    if len(steps) > 10:
        break

print(steps)