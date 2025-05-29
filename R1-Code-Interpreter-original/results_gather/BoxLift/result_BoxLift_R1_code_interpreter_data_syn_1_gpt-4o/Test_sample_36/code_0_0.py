boxes = [85, 162, 147, 83, 142, 96, 200, 172, 151, 77, 59, 39]
lifters = [53, 54, 76, 48, 97]

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
            # Box can be lifted, add to step
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