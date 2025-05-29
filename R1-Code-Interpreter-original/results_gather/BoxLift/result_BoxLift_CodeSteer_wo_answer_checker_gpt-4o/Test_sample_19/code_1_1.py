boxes = [55, 77, 274, 297, 70, 312, 108, 84, 340, 273, 157, 386, 95, 223, 311, 306, 287, 91, 220, 296, 320, 237, 97, 73]
lifters = [167, 138, 133, 87, 142, 167, 140]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
box_index = 0

while box_index < len(boxes) and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    for i in range(len(boxes)):
        if box_index >= len(boxes):
            break
        box_weight = boxes[box_index]
        lifter_indices = []
        total_capacity = 0
        for j in range(len(lifters)):
            if not used_lifters[j] and total_capacity < box_weight:
                lifter_indices.append(j)
                total_capacity += lifters[j]
                used_lifters[j] = True
            if total_capacity >= box_weight:
                break
        if total_capacity >= box_weight:
            step.append((box_weight, lifter_indices))
            box_index += 1
    steps.append(step)

# Check if all boxes are lifted within 7 steps
if box_index < len(boxes):
    print("It's not possible to lift all boxes in 7 or fewer steps.")
else:
    # Output the steps
    output = "<<<"
    for i, step in enumerate(steps):
        output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)