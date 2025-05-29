# Box weights and lifter capacities
boxes = [85, 162, 147, 83, 142, 96, 200, 172, 151, 77, 59, 39]
lifters = [53, 54, 76, 48, 97]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        remaining_lifters = lifters[:]
        for box in boxes[:]:
            # Try to find a single lifter or a combination of lifters to lift the box
            for i, lifter in enumerate(remaining_lifters):
                if lifter >= box:
                    step.append((box, [i]))
                    remaining_lifters.pop(i)
                    boxes.remove(box)
                    break
            else:
                # Try combinations of lifters
                for i in range(len(remaining_lifters)):
                    for j in range(i + 1, len(remaining_lifters)):
                        if remaining_lifters[i] + remaining_lifters[j] >= box:
                            step.append((box, [i, j]))
                            remaining_lifters.pop(j)
                            remaining_lifters.pop(i)
                            boxes.remove(box)
                            break
                    else:
                        continue
                    break
        steps.append(step)
        if len(steps) > 6:
            return "Not possible to lift all boxes in 6 steps"
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
print(steps)