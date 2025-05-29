boxes = [44, 178, 94, 176, 159, 199, 197, 133, 173, 186, 196, 80]
lifters = [61, 45, 66, 38, 52]

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
            # Assign the lifters to this box
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