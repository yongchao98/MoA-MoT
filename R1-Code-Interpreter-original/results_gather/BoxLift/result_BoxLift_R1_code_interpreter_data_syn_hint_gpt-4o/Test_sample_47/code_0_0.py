boxes = [280, 32, 168, 120, 112, 285, 46, 33, 148, 74, 295, 279, 207, 138, 153, 63]
lifters = [53, 92, 104, 86, 106, 68]

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
            if remaining_capacity <= 0:
                break
            if lifter <= remaining_capacity:
                used_lifters.append(i)
                remaining_capacity -= lifter
        if remaining_capacity <= 0:
            # Box can be lifted
            step.append((box, used_lifters))
            boxes.remove(box)
            # Remove used lifters
            for index in sorted(used_lifters, reverse=True):
                del remaining_lifters[index]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")