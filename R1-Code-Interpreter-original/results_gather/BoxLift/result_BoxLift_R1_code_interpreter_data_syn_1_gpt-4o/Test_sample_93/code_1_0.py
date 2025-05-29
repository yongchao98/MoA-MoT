boxes = [64, 70, 386, 351, 113, 77, 314, 333, 266, 399, 193, 44, 181, 200, 238, 175, 370, 118, 337, 134]
lifters = [140, 115, 159, 147, 129, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
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
        else:
            remaining_boxes.append(box)
    
    boxes = remaining_boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")