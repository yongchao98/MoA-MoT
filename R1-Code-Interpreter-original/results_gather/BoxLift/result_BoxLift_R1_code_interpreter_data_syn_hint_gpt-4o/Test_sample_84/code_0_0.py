boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
lifters = [164, 182, 169, 127, 118, 60]

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
            # Remove the used lifters from the remaining lifters
            for index in sorted(used_lifters, reverse=True):
                del remaining_lifters[index]
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")