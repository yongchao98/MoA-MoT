boxes = [45, 32, 84, 101, 161, 118, 30, 152, 159, 20, 110, 103]
lifters = [39, 32, 112, 102, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    # Try to find a combination of lifters to lift the box
    for i in range(len(lifters)):
        if used_lifters[i]:
            continue
        if lifters[i] >= box:
            return [i]
        # Try combinations
        total_capacity = lifters[i]
        lifter_indices = [i]
        for j in range(i + 1, len(lifters)):
            if used_lifters[j]:
                continue
            total_capacity += lifters[j]
            lifter_indices.append(j)
            if total_capacity >= box:
                return lifter_indices
    return None

def backtrack(step, used_boxes, steps):
    if all(used_boxes):
        return steps
    if step >= 5:
        return None

    used_lifters = [False] * len(lifters)
    step_lift = []

    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step_lift.append((box, lifter_indices))
            for idx in lifter_indices:
                used_lifters[idx] = True
            used_boxes[i] = True

    if step_lift:
        steps.append(step_lift)
        result = backtrack(step + 1, used_boxes, steps)
        if result is not None:
            return result
        steps.pop()
        for box, lifter_indices in step_lift:
            for idx in lifter_indices:
                used_lifters[idx] = False
            used_boxes[boxes.index(box)] = False

    return None

used_boxes = [False] * len(boxes)
steps = backtrack(0, used_boxes, [])

# Output the steps
if steps:
    output = ""
    for i, step in enumerate(steps):
        output += f"Step {i + 1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("<<<No solution found within 5 steps>>>")