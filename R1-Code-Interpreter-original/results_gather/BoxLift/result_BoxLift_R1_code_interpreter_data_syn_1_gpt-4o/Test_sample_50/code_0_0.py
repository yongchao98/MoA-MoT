# Box weights and lifter capacities
boxes = [92, 26, 84, 55, 72, 48, 33, 19]
lifters = [64, 48, 63]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Function to find the optimal lifting plan
def find_lifting_plan(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for i, box in enumerate(boxes[:]):
            # Try to lift the box with available lifters
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    boxes.remove(box)
                    break
            else:
                # Try to combine lifters if a single lifter can't lift the box
                for j in range(len(lifters)):
                    if used_lifters[j]:
                        continue
                    for k in range(j + 1, len(lifters)):
                        if used_lifters[k]:
                            continue
                        if lifters[j] + lifters[k] >= box:
                            step.append((box, [j, k]))
                            used_lifters[j] = True
                            used_lifters[k] = True
                            boxes.remove(box)
                            break
                    if used_lifters[j]:
                        break
        steps.append(step)
        if len(steps) > 4:
            return "Not possible to lift all boxes in 4 steps"
    return steps

# Find the lifting plan
lifting_plan = find_lifting_plan(boxes, lifters)

# Print the lifting plan
print(lifting_plan)