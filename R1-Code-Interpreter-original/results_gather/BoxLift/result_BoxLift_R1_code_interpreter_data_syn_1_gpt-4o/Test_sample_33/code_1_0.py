boxes = [88, 47, 209, 393, 59, 342, 259, 103, 304, 239, 141, 90, 256, 176, 297, 274, 154, 159, 262, 183]
lifters = [173, 188, 50, 119, 177, 86]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        # Try to find a single lifter or a combination of lifters to lift the box
        combination = []
        total_capacity = 0
        for i, lifter in enumerate(available_lifters):
            if not used_lifters[i] and total_capacity < box:
                combination.append(i)
                total_capacity += lifter
                used_lifters[i] = True
            if total_capacity >= box:
                break
        if total_capacity >= box:
            step.append((box, combination))
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")