boxes = [91, 207, 152, 47, 209, 54, 251, 176, 194, 221, 152, 141, 128, 159, 57, 184]
lifters = [149, 131, 113, 109, 124]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters to lift the box
    for j in range(1 << len(lifters)):
        total_capacity = 0
        lifter_indices = []
        for k in range(len(lifters)):
            if j & (1 << k) and not used[k]:
                total_capacity += lifters[k]
                lifter_indices.append(k)
        if total_capacity >= box:
            return lifter_indices
    return None

def find_steps(boxes, lifters, step, steps, max_steps, used):
    if not boxes:
        return True
    if step >= max_steps:
        return False

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            # Assign this combination to the current step
            steps[step].append((box, lifter_indices))
            # Mark lifters as used
            for idx in lifter_indices:
                used[idx] = True
            # Recurse to the next box
            if find_steps(boxes[:i] + boxes[i+1:], lifters, step + 1, steps, max_steps, used):
                return True
            # Backtrack
            steps[step].pop()
            for idx in lifter_indices:
                used[idx] = False
    return False

# Initialize steps and used lifters
steps = [[] for _ in range(6)]
used = [False] * len(lifters)

if find_steps(boxes, lifters, 0, steps, 6, used):
    output = []
    for i, step in enumerate(steps):
        if step:
            output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 6 steps.")