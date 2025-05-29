boxes = [87, 123, 302, 308, 370, 86, 182, 79, 135, 170, 198, 243, 45, 212, 183, 174, 335, 134, 58, 131, 124, 360, 192, 245]
lifters = [167, 156, 74, 177, 77, 60]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = boxes[:]
    
    for box in remaining_boxes:
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
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"
output += ">>>"
print(output)