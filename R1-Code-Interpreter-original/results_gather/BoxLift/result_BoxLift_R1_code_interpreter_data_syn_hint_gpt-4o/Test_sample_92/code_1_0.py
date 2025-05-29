boxes = [187, 170, 257, 238, 229, 202, 239, 125, 197, 243, 110, 50, 48, 168, 42, 292]
lifters = [115, 55, 142, 79, 63]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        remaining_capacity = box
        lifter_indices = []
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_capacity:
                lifter_indices.append(j)
                remaining_capacity -= lifter
                used_lifters[j] = True
                if remaining_capacity <= 0:
                    break
        
        if remaining_capacity <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)  # Remove the box from the original list
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")