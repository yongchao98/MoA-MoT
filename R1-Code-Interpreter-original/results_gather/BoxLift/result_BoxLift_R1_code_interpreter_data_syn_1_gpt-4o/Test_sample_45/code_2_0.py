boxes = [100, 22, 67, 79, 32, 31, 67, 37]
lifters = [58, 78, 68, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(3):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        lifter_indices = []
        remaining_weight = box
        
        # Try to find a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j]:
                lifter_indices.append(j)
                used_lifters[j] = True
                remaining_weight -= lifter
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step_lifting.append((box, lifter_indices))
            used_boxes[i] = True
    
    steps.append(step_lifting)
    
    if all(used_boxes):
        break

print(steps)