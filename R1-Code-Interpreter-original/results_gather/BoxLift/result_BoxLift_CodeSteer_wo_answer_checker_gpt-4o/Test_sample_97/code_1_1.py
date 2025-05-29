boxes = [184, 93, 275, 216, 137, 181, 31, 79, 56, 138, 81, 205, 108, 193, 230, 252]
lifters = [47, 157, 156, 45, 151, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity + lifter <= box:
            total_capacity += lifter
            lifter_indices.append(i)
            used[i] = True
            if total_capacity >= box:
                return lifter_indices
    return None

def assign_boxes(boxes, lifters, max_steps):
    steps = []
    for _ in range(max_steps):
        step = []
        used = [False] * len(lifters)
        for box in boxes[:]:
            lifter_indices = can_lift(box, lifters, used)
            if lifter_indices is not None:
                step.append((box, lifter_indices))
                boxes.remove(box)
        steps.append(step)
        if not boxes:
            break
    return steps

steps = assign_boxes(boxes, lifters, 6)

# Format the output
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")