boxes = [44, 70, 139, 184, 194, 67, 188, 33, 149, 194, 132, 117]
lifters = [110, 109, 30, 34]

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
        remaining_weight = box
        for i, lifter in enumerate(remaining_lifters):
            if lifter <= remaining_weight:
                used_lifters.append(i)
                remaining_weight -= lifter
                if remaining_weight <= 0:
                    break
        if remaining_weight <= 0:
            # If we can lift the box, add it to the step
            step.append((box, used_lifters))
            # Remove used lifters
            for index in sorted(used_lifters, reverse=True):
                del remaining_lifters[index]
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")