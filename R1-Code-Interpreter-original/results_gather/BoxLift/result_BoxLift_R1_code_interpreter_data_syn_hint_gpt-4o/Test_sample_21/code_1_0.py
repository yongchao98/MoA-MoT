boxes = [38, 97, 21, 92, 195, 187, 95, 140, 37, 85, 25, 165]
lifters = [96, 32, 88, 48]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):  # Iterate over a copy of the list
        # Try to find a single lifter or a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes.remove(box)
                break
        else:
            # Try to find a combination of lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                total_capacity = lifters[j]
                lifter_indices = [j]
                used_lifters[j] = True
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if total_capacity + lifters[k] >= box:
                        lifter_indices.append(k)
                        used_lifters[k] = True
                        step.append((box, lifter_indices))
                        boxes.remove(box)
                        break
                else:
                    # If no combination found, reset used lifters
                    for idx in lifter_indices:
                        used_lifters[idx] = False
                    continue
                break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")