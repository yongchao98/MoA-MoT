boxes = [174, 117, 74, 120, 98, 65, 112, 61, 46, 87, 176, 110]
lifters = [119, 104, 67, 98]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
max_steps = 5

while boxes and len(steps) < max_steps:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
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
            boxes.remove(box)  # Remove the box from the original list
    
    if step:
        steps.append(step)
    else:
        break  # If no boxes can be lifted in this step, break the loop

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")