boxes = [92, 26, 84, 55, 72, 48, 33, 19]
lifters = [64, 48, 63]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes.copy()

while remaining_boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to find a single lifter for the box
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and lifter >= box:
                step.append((box, [i]))
                used_lifters[i] = True
                remaining_boxes.remove(box)
                break
        else:
            # Try to find a combination of lifters for the box
            for i in range(len(lifters)):
                if used_lifters[i]:
                    continue
                for j in range(i + 1, len(lifters)):
                    if used_lifters[j]:
                        continue
                    if lifters[i] + lifters[j] >= box:
                        step.append((box, [i, j]))
                        used_lifters[i] = True
                        used_lifters[j] = True
                        remaining_boxes.remove(box)
                        break
                if used_lifters[i]:
                    break
    
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")