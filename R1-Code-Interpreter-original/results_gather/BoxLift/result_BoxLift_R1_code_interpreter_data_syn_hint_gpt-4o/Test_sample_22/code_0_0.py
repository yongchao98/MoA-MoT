boxes = [215, 113, 62, 255, 246, 148, 283, 149, 295, 250, 114, 48, 245, 145, 153, 50]
lifters = [82, 55, 159, 56, 88]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        used_lifters = []
        total_capacity = 0
        for i, lifter in enumerate(remaining_lifters):
            if total_capacity < box:
                used_lifters.append(i)
                total_capacity += lifter
            if total_capacity >= box:
                break
        if total_capacity >= box:
            # Assign the box to the lifters
            step.append((box, used_lifters))
            # Remove used lifters
            for i in sorted(used_lifters, reverse=True):
                del remaining_lifters[i]
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")