boxes = [45, 139, 291, 189, 171, 265, 139, 312, 289, 394, 183, 229, 321, 250, 226, 335, 59, 127, 298, 339, 256, 161, 288, 179]
lifters = [52, 99, 86, 56, 92, 134]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters to lift the box
    lifter_indices = []
    remaining_weight = box
    for j, lifter in enumerate(lifters):
        if not used[j] and lifter <= remaining_weight:
            lifter_indices.append(j)
            remaining_weight -= lifter
            used[j] = True
        if remaining_weight <= 0:
            return lifter_indices
    return None

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = boxes[:]
    
    for box in remaining_boxes:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")