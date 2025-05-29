boxes = [254, 205, 272, 360, 320, 328, 267, 223, 337, 238, 134, 235, 388, 160, 244, 378, 349, 262, 104, 396]
lifters = [132, 86, 81, 199, 93, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and remaining_weight > 0:
                if capacity <= remaining_weight:
                    lifter_indices.append(i)
                    remaining_weight -= capacity
                    used_lifters[i] = True
                if remaining_weight <= 0:
                    break
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")