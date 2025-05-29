boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
lifters = [164, 182, 169, 127, 118, 60]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    lifter_indices = []
    total_capacity = 0
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity + lifter <= box:
            total_capacity += lifter
            lifter_indices.append(i)
            used[i] = True
            if total_capacity >= box:
                return lifter_indices
    return None

def lift_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used = [False] * len(lifters)
        for box in boxes[:]:
            lifter_indices = can_lift(box, lifters, used)
            if lifter_indices is not None:
                step.append((box, lifter_indices))
                boxes.remove(box)
        steps.append(step)
    return steps

steps = lift_boxes(boxes, lifters)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")