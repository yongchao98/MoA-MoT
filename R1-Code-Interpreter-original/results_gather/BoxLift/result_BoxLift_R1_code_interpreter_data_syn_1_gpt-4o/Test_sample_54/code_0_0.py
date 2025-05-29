boxes = [87, 123, 302, 308, 370, 86, 182, 79, 135, 170, 198, 243, 45, 212, 183, 174, 335, 134, 58, 131, 124, 360, 192, 245]
lifters = [167, 156, 74, 177, 77, 60]

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
            # Remove the used lifters
            for i in sorted(used_lifters, reverse=True):
                del remaining_lifters[i]
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")