boxes = [156, 123, 139, 191, 91, 193, 100, 47, 36, 169, 90, 117]
lifters = [52, 55, 63, 44]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_indices = set()
    for box in boxes[:]:
        # Try to find a combination of lifters that can lift the box
        used_lifters = []
        total_capacity = 0
        for i, lifter in enumerate(remaining_lifters):
            if i not in used_indices and total_capacity < box:
                used_lifters.append(i)
                total_capacity += lifter
            if total_capacity >= box:
                break
        if total_capacity >= box:
            step.append((box, used_lifters))
            boxes.remove(box)
            used_indices.update(used_lifters)
    steps.append(step)

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")