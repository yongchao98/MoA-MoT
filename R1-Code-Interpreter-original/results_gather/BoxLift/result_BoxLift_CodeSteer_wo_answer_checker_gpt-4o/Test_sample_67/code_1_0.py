def box_lift(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        
        for i, box in enumerate(boxes):
            if not boxes[i]:
                continue
            for j, lifter in enumerate(lifters):
                if used_lifters[j]:
                    continue
                if lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    boxes[i] = 0
                    break
            else:
                # Try to combine lifters
                total_capacity = 0
                lifter_indices = []
                for j, lifter in enumerate(lifters):
                    if used_lifters[j]:
                        continue
                    if total_capacity + lifter >= box:
                        lifter_indices.append(j)
                        step.append((box, lifter_indices))
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        boxes[i] = 0
                        break
                    else:
                        total_capacity += lifter
                        lifter_indices.append(j)
        
        steps.append(step)
        boxes = [b for b in boxes if b > 0]
        if len(steps) > 6:
            break
    
    return steps

boxes = [93, 50, 31, 217, 183, 34, 43, 268, 281, 93, 145, 74, 278, 272, 86, 81]
lifters = [97, 53, 143, 47, 94, 132]

steps = box_lift(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)