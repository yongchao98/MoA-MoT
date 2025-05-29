boxes = [247, 236, 255, 64, 58, 144, 218, 200, 320, 87, 256, 111, 379, 318, 367, 88, 136, 125, 172, 93]
lifters = [63, 51, 130, 110, 117, 132]

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
        remaining_capacity = box
        for i, lifter in enumerate(remaining_lifters):
            if lifter <= remaining_capacity:
                used_lifters.append(i)
                remaining_capacity -= lifter
                if remaining_capacity <= 0:
                    break
        if remaining_capacity <= 0:
            # If we found a valid combination, add it to the step
            step.append((box, used_lifters))
            # Remove the used lifters and the box
            for i in sorted(used_lifters, reverse=True):
                del remaining_lifters[i]
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")