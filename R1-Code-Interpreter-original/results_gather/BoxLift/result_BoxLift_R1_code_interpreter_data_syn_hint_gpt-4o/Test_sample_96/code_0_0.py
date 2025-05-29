boxes = [169, 106, 153, 20, 138, 80, 94, 190, 142, 196, 186, 87]
lifters = [55, 85, 43, 100, 47]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_lifters = [False] * len(lifters)

def can_lift(box_weight, lifters, used_lifters):
    # Try to find a combination of lifters that can lift the box
    current_weight = 0
    lifter_indices = []
    for i, capacity in enumerate(lifters):
        if not used_lifters[i] and current_weight < box_weight:
            current_weight += capacity
            lifter_indices.append(i)
            used_lifters[i] = True
        if current_weight >= box_weight:
            return lifter_indices
    return None

for _ in range(7):  # Maximum 7 steps
    step = []
    for box_weight in boxes:
        lifter_indices = can_lift(box_weight, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box_weight, lifter_indices))
            for idx in lifter_indices:
                used_lifters[idx] = True
    if step:
        steps.append(step)
    if len(steps) >= 7 or all(used_lifters):
        break
    used_lifters = [False] * len(lifters)  # Reset lifters for the next step

print(steps)