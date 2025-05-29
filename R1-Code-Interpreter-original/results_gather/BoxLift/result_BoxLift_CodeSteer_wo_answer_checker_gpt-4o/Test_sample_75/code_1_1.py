boxes = [76, 97, 14, 30, 36, 83, 45, 46]
lifters = [69, 59, 68, 50]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

for step in range(3):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = []
        remaining_weight = box
        
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and capacity <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= capacity
                used_lifters[i] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step_lift.append((box, lifter_indices))
            boxes.remove(box)
    
    if step_lift:
        steps.append(f"Step {step + 1}: {step_lift}")
    if not boxes:
        break

result = "\n".join(steps)
print(f"<<<{result}>>>")