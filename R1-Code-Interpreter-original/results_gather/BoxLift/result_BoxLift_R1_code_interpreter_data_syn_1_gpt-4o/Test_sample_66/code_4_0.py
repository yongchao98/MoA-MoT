boxes = [181, 52, 82, 130, 92, 195, 101, 106, 189, 42, 196, 46]
lifters = [45, 30, 52, 104, 30]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 8:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Use a copy of the list to iterate
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
    
    if step:
        steps.append(step)
    else:
        break  # If no boxes can be lifted in this step, break the loop

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")