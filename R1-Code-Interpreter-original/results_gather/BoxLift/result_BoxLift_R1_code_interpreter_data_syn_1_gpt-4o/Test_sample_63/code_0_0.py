boxes = [49, 186, 267, 243, 352, 74, 160, 115, 138, 301, 250, 145, 294, 232, 144, 293, 287, 358, 267, 266]
lifters = [60, 92, 94, 68, 72, 79]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        used_lifters = []
        remaining_capacity = box
        for i, lifter in enumerate(available_lifters):
            if remaining_capacity <= 0:
                break
            if lifter <= remaining_capacity:
                used_lifters.append(i)
                remaining_capacity -= lifter
        if remaining_capacity <= 0:
            # Box can be lifted
            step.append((box, used_lifters))
            # Remove used lifters from available lifters
            for index in sorted(used_lifters, reverse=True):
                del available_lifters[index]
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")