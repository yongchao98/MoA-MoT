boxes = [76, 97, 14, 30, 36, 83, 45, 46]
lifters = [69, 59, 68, 50]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters that can lift the box
    lifter_indices = []
    remaining_weight = box
    for i, capacity in enumerate(lifters):
        if not used[i] and capacity <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= capacity
            used[i] = True
            if remaining_weight <= 0:
                return lifter_indices
    # If we can't lift the box, reset the used lifters
    for i in lifter_indices:
        used[i] = False
    return None

steps = []

for step in range(3):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step_lift.append((box, lifter_indices))
            boxes.remove(box)
    
    if step_lift:
        steps.append(f"Step {step + 1}: {step_lift}")
    if not boxes:
        break

result = "\n".join(steps)
print(f"<<<{result}>>>")