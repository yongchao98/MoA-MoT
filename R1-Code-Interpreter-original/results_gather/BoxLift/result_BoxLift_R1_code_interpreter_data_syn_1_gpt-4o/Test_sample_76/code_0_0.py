boxes = [94, 177, 88, 131, 87, 113, 148, 86, 141, 83, 86, 73]
lifters = [55, 92, 35, 48, 96]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

for step in range(6):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        total_capacity = 0
        
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and total_capacity < box:
                lifter_indices.append(i)
                total_capacity += capacity
                used_lifters[i] = True
            
            if total_capacity >= box:
                break
        
        if total_capacity >= box:
            step_lift.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    steps.append(step_lift)
    
    if not remaining_boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")