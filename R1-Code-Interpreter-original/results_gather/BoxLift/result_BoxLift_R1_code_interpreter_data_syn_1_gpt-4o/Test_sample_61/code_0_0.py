boxes = [128, 144, 42, 49, 132, 118, 32, 48, 21, 94, 85, 51]
lifters = [85, 90, 40, 102]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes.copy()

for step in range(5):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to lift the box with available lifters
        lifter_indices = []
        remaining_weight = box
        
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and remaining_weight > 0:
                if capacity >= remaining_weight:
                    lifter_indices.append(i)
                    used_lifters[i] = True
                    remaining_weight = 0
                    break
                else:
                    lifter_indices.append(i)
                    used_lifters[i] = True
                    remaining_weight -= capacity
        
        if remaining_weight == 0:
            step_lift.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    steps.append(step_lift)
    
    if not remaining_boxes:
        break

print(steps)