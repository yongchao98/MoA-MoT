boxes = [16, 56, 15, 82, 85, 12, 12, 32]
lifters = [63, 48, 72, 68]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(3):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if used_lifters[j]:
                continue
            if lifter >= box:
                # Single lifter can lift the box
                step_lift.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to find a combination of lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        # Two lifters can lift the box
                        step_lift.append((box, [j, k]))
                        used_boxes[i] = True
                        used_lifters[j] = True
                        used_lifters[k] = True
                        break
                if used_boxes[i]:
                    break
    
    steps.append(step_lift)
    if all(used_boxes):
        break

# Format the output
output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"
output += ">>>"

print(output)