boxes = [116, 246, 296, 369, 78, 275, 77, 383, 71, 155, 352, 94, 63, 168, 350, 79, 59, 252, 88, 278, 188, 383, 240, 308]
lifters = [134, 132, 150, 137, 62, 132]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_capacity = lifters[:]
    
    for i, box in enumerate(boxes):
        if box is None:
            continue
        lifter_indices = []
        total_capacity = 0
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
                if total_capacity >= box:
                    break
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes[i] = None
    
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")