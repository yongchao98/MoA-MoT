boxes = [202, 297, 200, 241, 156, 177, 70, 78, 133, 176, 75, 70, 167, 240, 299, 91]
lifters = [56, 160, 118, 138, 160, 150]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(5):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a combination of lifters to lift this box
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and remaining_weight > 0:
                lifter_indices.append(j)
                used_lifters[j] = True
                remaining_weight -= lifter
        
        if remaining_weight <= 0:
            step_lifting.append((box, lifter_indices))
            used_boxes[i] = True
    
    if not step_lifting:
        break
    
    steps.append(step_lifting)

# Check if all boxes are lifted
all_lifted = all(used_boxes)

if all_lifted:
    print(steps)
else:
    print("Not all boxes could be lifted within 5 steps.")