boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
lifters = [122, 116, 51, 165, 95, 150, 179]

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

# Check if all boxes are lifted
if box_index < len(boxes):
    print("Unable to lift all boxes within 7 steps.")
else:
    # Print the steps
    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)