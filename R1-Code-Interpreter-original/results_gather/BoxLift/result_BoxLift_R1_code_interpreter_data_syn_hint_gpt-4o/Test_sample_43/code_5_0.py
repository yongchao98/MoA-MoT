boxes = [247, 236, 255, 64, 58, 144, 218, 200, 320, 87, 256, 111, 379, 318, 367, 88, 136, 125, 172, 93]
lifters = [63, 51, 130, 110, 117, 132]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 9:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        used_indices = []
        remaining_capacity = box
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and lifter <= remaining_capacity:
                used_indices.append(i)
                remaining_capacity -= lifter
                used_lifters[i] = True
                if remaining_capacity <= 0:
                    break
        if remaining_capacity <= 0:
            step.append((box, used_indices))
            boxes.remove(box)
    if step:
        steps.append(step)
    else:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")