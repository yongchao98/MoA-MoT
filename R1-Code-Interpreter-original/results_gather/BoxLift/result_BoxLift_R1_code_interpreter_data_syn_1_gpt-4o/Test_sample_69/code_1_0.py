boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

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
    
    steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")