boxes = [14, 51, 14, 69, 77, 20, 52, 96]
lifters = [69, 66, 74]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(4):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step_lift.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to combine lifters to lift the box
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    if total_capacity >= box:
                        step_lift.append((box, lifter_indices))
                        used_boxes[i] = True
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        break
    
    steps.append(step_lift)
    if all(used_boxes):
        break

output = "\n".join(f"Step {i+1}: {step}" for i, step in enumerate(steps))
print(f"<<<{output}>>>")