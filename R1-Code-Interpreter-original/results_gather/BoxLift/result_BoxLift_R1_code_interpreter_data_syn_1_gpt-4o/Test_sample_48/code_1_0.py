boxes = [295, 145, 158, 66, 383, 291, 156, 366, 271, 121, 289, 388, 276, 247, 68, 248, 167, 252, 126, 159]
lifters = [144, 152, 155, 109, 115, 152, 162]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(6):
    step_lifts = []
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
                remaining_weight -= lifter
                used_lifters[j] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step_lifts.append((box, lifter_indices))
            used_boxes[i] = True
    
    steps.append(step_lifts)
    
    # Check if all boxes are lifted
    if all(used_boxes):
        break

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")