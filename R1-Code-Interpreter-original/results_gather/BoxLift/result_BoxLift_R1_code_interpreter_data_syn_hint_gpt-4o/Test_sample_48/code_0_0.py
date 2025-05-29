boxes = [295, 145, 158, 66, 383, 291, 156, 366, 271, 121, 289, 388, 276, 247, 68, 248, 167, 252, 126, 159]
lifters = [144, 152, 155, 109, 115, 152, 162]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

while remaining_boxes and len(steps) < 6:
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
    print(f"Step {i + 1}: {step}")