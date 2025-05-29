boxes = [79, 48, 16, 95, 67, 41, 62, 22]
lifters = [52, 41, 41, 67]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(3):
    step_lifts = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lifts.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to combine lifters for the box
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    if total_capacity >= box:
                        step_lifts.append((box, lifter_indices))
                        used_boxes[i] = True
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        break
    
    steps.append(step_lifts)
    if all(used_boxes):
        break

# Format the output
output = []
for i, step in enumerate(steps):
    step_str = f"Step {i + 1}: {step}"
    output.append(step_str)

result = "\n".join(output)
print(f"<<<{result}>>>")