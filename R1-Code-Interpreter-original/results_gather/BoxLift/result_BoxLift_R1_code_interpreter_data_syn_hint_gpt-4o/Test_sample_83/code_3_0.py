# Box weights and lifter capacities
boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes and len(steps) < 7:
        step = []
        used_lifters = [False] * len(lifters)
        for box in boxes[:]:
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step.append((box, lifter_indices))
                        boxes.remove(box)
                        break
        steps.append(step)
    if boxes:
        return "Cannot lift all boxes in 7 steps"
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
print(steps)