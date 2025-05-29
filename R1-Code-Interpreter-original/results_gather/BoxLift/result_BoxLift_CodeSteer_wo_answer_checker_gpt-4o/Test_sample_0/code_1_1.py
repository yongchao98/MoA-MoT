boxes = [18, 93, 23, 70, 72, 22, 96, 96]
lifters = [60, 46, 61]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_capacity = lifters[:]
    
    for i, box in enumerate(boxes):
        if box is None:
            continue
        
        # Try to find a single lifter or a combination of lifters to lift the box
        lifter_indices = []
        total_capacity = 0
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
                if total_capacity >= box:
                    break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes[i] = None
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"
print(output)