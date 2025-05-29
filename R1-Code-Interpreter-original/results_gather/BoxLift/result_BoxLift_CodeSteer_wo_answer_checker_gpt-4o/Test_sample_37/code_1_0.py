boxes = [97, 164, 187, 145, 255, 257, 240, 55, 261, 185, 252, 257, 206, 298, 121, 202]
lifters = [131, 106, 144, 121, 109, 52]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    used_lifters = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            used_lifters.append(i)
            total_capacity += lifter
        if total_capacity >= box:
            return used_lifters
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 7:
        return False

    for i, box in enumerate(boxes):
        used_lifters = can_lift(box, lifters, used)
        if used_lifters is not None:
            # Mark lifters as used
            for index in used_lifters:
                used[index] = True
            # Try to lift the next box
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step, steps, used):
                steps[step].append((box, used_lifters))
                return True
            # Backtrack
            for index in used_lifters:
                used[index] = False

    # Move to the next step
    if backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters)):
        return True

    return False

steps = [[] for _ in range(7)]
used = [False] * len(lifters)
if backtrack(boxes, lifters, 0, steps, used):
    output = []
    for idx, step in enumerate(steps):
        if step:
            output.append(f"Step {idx + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 7 steps.")