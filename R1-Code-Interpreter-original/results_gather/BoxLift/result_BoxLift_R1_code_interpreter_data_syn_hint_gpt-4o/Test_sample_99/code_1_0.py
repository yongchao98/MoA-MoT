boxes = [289, 375, 107, 145, 257, 48, 141, 83, 136, 368, 59, 133, 186, 266, 353, 73, 66, 210, 247, 79, 342, 318, 337, 162]
lifters = [75, 189, 162, 181, 118, 194, 137]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

steps = []

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        lifter_indices = []
        remaining_weight = box
        
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and lifter <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= lifter
                used_lifters[i] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")