boxes = [82, 61, 91, 29, 20, 43, 57, 79]
lifters = [49, 80, 79]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to lift the box with available lifters
    for i in range(len(lifters)):
        if not used[i] and lifters[i] >= box:
            return [i]
    # Try combinations of two lifters
    for i in range(len(lifters)):
        if used[i]:
            continue
        for j in range(i + 1, len(lifters)):
            if used[j]:
                continue
            if lifters[i] + lifters[j] >= box:
                return [i, j]
    return None

def backtrack(steps, current_step, boxes, lifters, max_steps):
    if not boxes:
        return steps if len(steps) <= max_steps else None
    if current_step >= max_steps:
        return None

    used = [False] * len(lifters)
    step = []
    remaining_boxes = boxes[:]
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            for idx in lifter_indices:
                used[idx] = True
            remaining_boxes.remove(box)

    if step:
        steps.append(step)
        result = backtrack(steps, current_step + 1, remaining_boxes, lifters, max_steps)
        if result:
            return result
        steps.pop()

    return None

# Start backtracking
solution = backtrack([], 0, boxes, lifters, 4)

# Output the solution
if solution:
    output = ""
    for i, step in enumerate(solution):
        output += f"Step {i + 1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("No solution found within 4 steps.")