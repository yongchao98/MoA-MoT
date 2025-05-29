# Box weights and lifter capacities
boxes = [24, 78, 95, 23, 57, 36, 40, 26]
lifters = [56, 51, 66, 65]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Function to find the minimum number of steps to lift all boxes
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
                # Try to combine lifters to lift the box
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
        if len(steps) > 3:
            return "Not possible to lift all boxes in 3 steps"
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
print(steps)