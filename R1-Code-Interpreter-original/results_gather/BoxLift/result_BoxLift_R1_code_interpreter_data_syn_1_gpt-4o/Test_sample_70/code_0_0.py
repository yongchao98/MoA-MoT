# Box weights and lifter capacities
boxes = [82, 61, 91, 29, 20, 43, 57, 79]
lifters = [49, 80, 79]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Function to find the optimal assignment of lifters to boxes
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
    return steps

# Get the steps for lifting all boxes
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")