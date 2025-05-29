boxes = [88, 47, 209, 393, 59, 342, 259, 103, 304, 239, 141, 90, 256, 176, 297, 274, 154, 159, 262, 183]
lifters = [173, 188, 50, 119, 177, 86]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = []
        total_capacity = 0
        
        # Try to lift the box with available lifters
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and total_capacity + capacity <= box:
                total_capacity += capacity
                lifter_indices.append(i)
                used_lifters[i] = True
                if total_capacity >= box:
                    break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)
    if len(steps) > 8:
        break

# Print the steps
output = []
for i, step in enumerate(steps):
    if step:
        output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")