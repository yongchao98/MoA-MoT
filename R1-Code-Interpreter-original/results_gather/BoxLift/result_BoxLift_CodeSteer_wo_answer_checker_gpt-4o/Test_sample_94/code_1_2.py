from functools import lru_cache

boxes = [187, 344, 112, 66, 384, 247, 184, 298, 131, 51, 62, 255, 120, 357, 399, 287, 231, 161, 336, 256, 328, 239, 365, 245]
lifters = [165, 85, 52, 82, 114, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

@lru_cache(None)
def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            total_capacity += lifter
            lifter_indices.append(i)
        if total_capacity >= box:
            return tuple(lifter_indices)
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 11:
        return False

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, tuple(lifters), tuple(used))
        if lifter_indices is not None:
            # Mark lifters as used
            for index in lifter_indices:
                used[index] = True
            # Add to current step
            steps[step].append((box, list(lifter_indices)))
            # Recurse to next box
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step, steps, used):
                return True
            # Backtrack
            for index in lifter_indices:
                used[index] = False
            steps[step].pop()

    # Move to next step
    if steps[step]:
        return backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters))
    return False

steps = [[] for _ in range(11)]
used = [False] * len(lifters)

if backtrack(tuple(boxes), tuple(lifters), 0, steps, used):
    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within 11 steps.")