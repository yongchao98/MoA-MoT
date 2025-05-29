boxes = [202, 297, 200, 241, 156, 177, 70, 78, 133, 176, 75, 70, 167, 240, 299, 91]
lifters = [56, 160, 118, 138, 160, 150]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        total_capacity = 0
        lifter_indices = []
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
                if total_capacity >= box:
                    break
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    steps.append(step)

# Format the output
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")