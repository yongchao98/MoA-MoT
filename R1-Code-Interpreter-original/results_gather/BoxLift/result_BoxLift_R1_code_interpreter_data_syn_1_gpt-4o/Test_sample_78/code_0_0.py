boxes = [260, 319, 78, 313, 164, 297, 367, 392, 175, 245, 200, 45, 386, 124, 169, 97, 359, 273, 175, 238]
lifters = [57, 82, 155, 87, 174, 192, 149]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    for box in boxes[:]:
        # Try to find a combination of lifters that can lift the box
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
            # Remove the used lifters from the remaining lifters
            for index in sorted(used_lifters, reverse=True):
                del remaining_lifters[index]
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)
    if len(steps) > 7:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")