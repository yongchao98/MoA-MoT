boxes = [119, 110, 359, 148, 237, 71, 309, 213, 117, 394, 79, 353, 230, 77, 242, 239, 66, 259, 212, 358, 194, 397, 274, 101]
lifters = [60, 103, 156, 108, 179, 110]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = []
        remaining_weight = box
        
        # Try to find a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                remaining_weight -= lifter
                used_lifters[j] = True
            if remaining_weight <= 0:
                break
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)  # Remove the lifted box
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")