boxes = [55, 244, 173, 293, 90, 126, 340, 250, 66, 143, 103, 244, 76, 166, 130, 216, 54, 196, 245, 307]
lifters = [142, 178, 196, 52, 101, 144, 50]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes[:]):
        total_capacity = 0
        lifter_indices = []
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity < box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
            if total_capacity >= box:
                break
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")