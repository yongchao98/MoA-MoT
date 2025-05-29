boxes = [58, 45, 94, 26, 49, 153, 61, 23, 170, 143, 75, 166]
lifters = [57, 61, 104, 98, 70]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = set()

for step in range(4):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if i in used_boxes:
            continue
        
        # Try to find a single lifter or a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                used_lifters[j] = True
                remaining_weight -= lifter
            
            if remaining_weight <= 0:
                break
        
        if remaining_weight <= 0:
            step_lifting.append((box, lifter_indices))
            used_boxes.add(i)
    
    steps.append(step_lifting)
    
    if len(used_boxes) == len(boxes):
        break

print(steps)