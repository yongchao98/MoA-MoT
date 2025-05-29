boxes = [17, 69, 16, 83, 95, 91, 68, 25]
lifters = [42, 68, 50, 48]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for _ in range(4):  # We have a maximum of 4 steps
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a combination of lifters to lift this box
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
            used_boxes[i] = True
    
    if step:
        steps.append(step)
    else:
        break

print(steps)