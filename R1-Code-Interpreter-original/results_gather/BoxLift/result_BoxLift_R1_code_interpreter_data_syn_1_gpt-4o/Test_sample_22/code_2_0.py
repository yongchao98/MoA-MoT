boxes = [215, 113, 62, 255, 246, 148, 283, 149, 295, 250, 114, 48, 245, 145, 153, 50]
lifters = [82, 55, 159, 56, 88]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        total_capacity = 0
        for i, lifter in enumerate(remaining_lifters):
            if not used_lifters[i] and total_capacity < box:
                lifter_indices.append(i)
                total_capacity += lifter
                used_lifters[i] = True
            if total_capacity >= box:
                break
        if total_capacity >= box:
            step.append((box, lifter_indices))
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")