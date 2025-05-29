boxes = [280, 109, 195, 129, 54, 216, 157, 193, 44, 233, 171, 218, 204, 120, 131, 151]
lifters = [49, 131, 66, 124, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        # Try to find a combination of lifters that can lift the box
        used_lifters = []
        remaining_weight = box
        for i, lifter in enumerate(available_lifters):
            if remaining_weight <= 0:
                break
            if lifter <= remaining_weight:
                used_lifters.append(i)
                remaining_weight -= lifter
        if remaining_weight <= 0:
            # If we found a valid combination, assign the lifters to the box
            step.append((box, used_lifters))
            boxes.remove(box)
            # Remove used lifters from available lifters
            for index in sorted(used_lifters, reverse=True):
                del available_lifters[index]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")