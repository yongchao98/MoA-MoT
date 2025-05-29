boxes = [181, 52, 82, 130, 92, 195, 101, 106, 189, 42, 196, 46]
lifters = [45, 30, 52, 104, 30]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    remaining_weight = box
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            used[i] = True
            if remaining_weight <= 0:
                return lifter_indices
    return None

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Use a copy of the list to iterate
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)  # Remove the box once it's lifted
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")