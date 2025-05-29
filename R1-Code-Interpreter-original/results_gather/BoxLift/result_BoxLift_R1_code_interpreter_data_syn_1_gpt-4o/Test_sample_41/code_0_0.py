boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
lifters = [95, 111, 115, 35]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(5):
    step_lifts = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
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
            step_lifts.append((box, lifter_indices))
            used_boxes[i] = True
    
    if step_lifts:
        steps.append(step_lifts)
    else:
        break

print(steps)