boxes = [75, 298, 70, 200, 278, 35, 278, 56, 270, 212, 161, 98, 244, 292, 177, 145]
lifters = [105, 102, 125, 93, 123]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        total_capacity = 0
        lifter_indices = []
        for j, capacity in enumerate(lifters):
            if not used_lifters[j] and total_capacity + capacity <= box:
                total_capacity += capacity
                lifter_indices.append(j)
                used_lifters[j] = True
            if total_capacity >= box:
                break
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")