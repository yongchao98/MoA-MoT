boxes = [57, 112, 121, 49, 166, 143, 169, 43, 90, 66, 23, 143]
lifters = [103, 49, 95, 97, 98]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_lifters = [False] * len(lifters)

def can_lift(box_weight, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box_weight

def find_lifters_for_box(box_weight):
    lifter_indices = []
    for i, capacity in enumerate(lifters):
        if not used_lifters[i]:
            lifter_indices.append(i)
            if can_lift(box_weight, lifter_indices):
                for idx in lifter_indices:
                    used_lifters[idx] = True
                return lifter_indices
    return None

for step in range(4):
    step_lift = []
    used_lifters = [False] * len(lifters)
    for box in boxes:
        if not any(box in s for s in steps):
            lifter_indices = find_lifters_for_box(box)
            if lifter_indices is not None:
                step_lift.append((box, lifter_indices))
    if step_lift:
        steps.append(step_lift)
    if len(steps) == 4 or all(any(box in s for s in steps) for box in boxes):
        break

print(steps)