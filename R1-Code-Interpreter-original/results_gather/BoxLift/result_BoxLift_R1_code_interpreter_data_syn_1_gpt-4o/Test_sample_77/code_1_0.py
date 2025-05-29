boxes = [253, 162, 171, 57, 68, 142, 196, 83, 94, 287, 256, 236, 153, 79, 111, 42]
lifters = [61, 62, 122, 75, 143, 131]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):  # Use a copy of the list to iterate
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
            boxes.remove(box)  # Remove the box once it's lifted
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")