boxes = [71, 55, 11, 10, 71, 91, 14, 69]
lifters = [46, 54, 47]

# Sort boxes in descending order to try to lift the heaviest first
boxes.sort(reverse=True)

# Sort lifters in descending order to use the strongest lifters first
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(5):  # We have a maximum of 5 steps
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter who can lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lifting.append((box, [j]))
                used_lifters[j] = True
                used_boxes[i] = True
                break
        else:
            # Try to combine lifters to lift the box
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        step_lifting.append((box, [j, k]))
                        used_lifters[j] = True
                        used_lifters[k] = True
                        used_boxes[i] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(step_lifting)
    
    # Check if all boxes are lifted
    if all(used_boxes):
        break

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")