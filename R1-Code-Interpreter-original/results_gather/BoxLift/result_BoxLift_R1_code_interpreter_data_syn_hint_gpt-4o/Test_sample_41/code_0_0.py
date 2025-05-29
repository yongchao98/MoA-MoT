boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
lifters = [95, 111, 115, 35]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

while remaining_boxes:
    step = []
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
            step.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")