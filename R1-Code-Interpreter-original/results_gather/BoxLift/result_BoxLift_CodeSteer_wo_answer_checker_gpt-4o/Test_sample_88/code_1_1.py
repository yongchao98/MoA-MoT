boxes = [88, 313, 344, 387, 392, 370, 99, 389, 340, 391, 131, 147, 334, 251, 157, 144, 290, 167, 141, 79]
lifters = [148, 193, 54, 146, 106, 53, 88]

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
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
                if total_capacity == box:
                    break
        if total_capacity == box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")