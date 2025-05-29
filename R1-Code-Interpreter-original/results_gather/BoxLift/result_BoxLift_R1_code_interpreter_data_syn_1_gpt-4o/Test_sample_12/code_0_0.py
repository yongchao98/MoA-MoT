boxes = [50, 219, 120, 245, 204, 172, 274, 225, 94, 270, 158, 253, 174, 40, 252, 124]
lifters = [54, 61, 106, 60, 63]

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
            # Assign this box to the lifters
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