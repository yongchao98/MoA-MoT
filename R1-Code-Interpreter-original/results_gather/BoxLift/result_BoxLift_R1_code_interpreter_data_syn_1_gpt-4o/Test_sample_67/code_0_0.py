boxes = [93, 50, 31, 217, 183, 34, 43, 268, 281, 93, 145, 74, 278, 272, 86, 81]
lifters = [97, 53, 143, 47, 94, 132]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
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
                        step.append((box, lifter_indices))
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        break
            else:
                continue
        boxes[i] = None
    
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")