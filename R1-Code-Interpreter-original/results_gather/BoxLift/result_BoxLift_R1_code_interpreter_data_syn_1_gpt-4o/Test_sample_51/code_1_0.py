boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
lifters = [93, 71, 82, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes.copy()

for step in range(5):
    step_lifts = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to lift the box with available lifters
        lifter_indices = []
        remaining_capacity = box
        
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and remaining_capacity > 0:
                lifter_indices.append(i)
                remaining_capacity -= capacity
                used_lifters[i] = True
                if remaining_capacity <= 0:
                    break
        
        if remaining_capacity <= 0:
            step_lifts.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    steps.append(step_lifts)
    
    if not remaining_boxes:
        break

print(steps)