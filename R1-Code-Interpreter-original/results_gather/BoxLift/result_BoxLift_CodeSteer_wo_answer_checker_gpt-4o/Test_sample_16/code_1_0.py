def box_lift(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        
        for i, box in enumerate(boxes[:]):
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter >= box:
                    step.append((box, [j]))
                    used_lifters[j] = True
                    boxes.remove(box)
                    break
            else:
                # Try to combine lifters if a single one can't lift the box
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
                        total_capacity += lifters[k]
                        lifter_indices.append(k)
                    if box not in boxes:
                        break
        
        steps.append(step)
        if len(steps) > 5:
            return "Not possible to lift all boxes in 5 steps or less."
    
    return steps

boxes = [166, 194, 27, 144, 198, 172, 81, 117, 154, 183, 27, 184]
lifters = [97, 90, 80, 111, 76]

result = box_lift(boxes, lifters)
output = "\n".join([f"Step {i+1}: {step}" for i, step in enumerate(result)])
print(f"<<<{output}>>>")