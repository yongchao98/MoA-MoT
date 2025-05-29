boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
lifters = [78, 98, 51, 75, 43]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    # Try to find a combination of lifters to lift the box
    for i in range(len(lifters)):
        if used_lifters[i]:
            continue
        total_capacity = lifters[i]
        lifter_indices = [i]
        used_lifters[i] = True
        for j in range(i + 1, len(lifters)):
            if used_lifters[j]:
                continue
            if total_capacity + lifters[j] >= box:
                lifter_indices.append(j)
                used_lifters[j] = True
                return lifter_indices
            total_capacity += lifters[j]
            lifter_indices.append(j)
            used_lifters[j] = True
        for index in lifter_indices:
            used_lifters[index] = False
    return None

def backtrack(boxes, lifters, step, max_steps, used_boxes, steps):
    if all(used_boxes):
        return steps
    if step >= max_steps:
        return None
    
    current_step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            used_boxes[i] = True
    
    if current_step:
        steps.append(current_step)
        result = backtrack(boxes, lifters, step + 1, max_steps, used_boxes, steps)
        if result is not None:
            return result
        steps.pop()
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        used_boxes[i] = False
    
    return None

steps = backtrack(boxes, lifters, 0, 5, [False] * len(boxes), [])
print(steps)