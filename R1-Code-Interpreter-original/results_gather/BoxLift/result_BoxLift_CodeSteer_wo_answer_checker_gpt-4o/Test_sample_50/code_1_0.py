def box_lift():
    boxes = [92, 26, 84, 55, 72, 48, 33, 19]
    lifters = [64, 48, 63]

    # Sort boxes and lifters in descending order
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
                    for k in range(j + 1, len(lifters)):
                        if used_lifters[k]:
                            continue
                        if lifters[j] + lifters[k] >= box:
                            step.append((box, [j, k]))
                            used_lifters[j] = True
                            used_lifters[k] = True
                            boxes.remove(box)
                            break
                    if used_lifters[j]:
                        break

        steps.append(step)
        if len(steps) > 4:
            break

    # Print the steps
    result = ""
    for i, step in enumerate(steps):
        result += f"Step {i + 1}: {step}\n"
    print(f"<<<{result.strip()}>>>")

box_lift()