# Box weights and lifter capacities
boxes = [41, 69, 80, 19, 54, 85, 64, 35]
lifters = [45, 45, 62, 74]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        available_lifters = lifters[:]
        used_lifters = [False] * len(lifters)
        
        for i, box in enumerate(boxes[:]):
            for j, lifter in enumerate(available_lifters):
                if not used_lifters[j] and lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    boxes.remove(box)
                    break
            else:
                # Try to combine lifters if a single one can't lift the box
                for j in range(len(available_lifters)):
                    if used_lifters[j]:
                        continue
                    for k in range(j + 1, len(available_lifters)):
                        if used_lifters[k]:
                            continue
                        if available_lifters[j] + available_lifters[k] >= box:
                            step.append((box, [j, k]))
                            used_lifters[j] = True
                            used_lifters[k] = True
                            boxes.remove(box)
                            break
                    if used_lifters[j]:
                        break
        
        steps.append(step)
        if len(steps) > 4:
            break
    
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")