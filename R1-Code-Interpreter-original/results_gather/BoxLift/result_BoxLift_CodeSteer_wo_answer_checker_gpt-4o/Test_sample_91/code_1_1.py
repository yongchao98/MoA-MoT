boxes = [159, 180, 64, 184, 110, 43, 128, 47, 47, 138, 62, 97]
lifters = [81, 98, 62, 46]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    remaining_weight = box
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
        if remaining_weight <= 0:
            return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 7:
        return False

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            # Mark lifters as used
            for idx in lifter_indices:
                used[idx] = True
            # Add to current step
            steps[step].append((box, lifter_indices))
            # Recurse to next box
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step, steps, used):
                return True
            # Backtrack
            for idx in lifter_indices:
                used[idx] = False
            steps[step].pop()

    # Move to next step
    if step < 6:
        return backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters))
    return False

steps = [[] for _ in range(7)]
used = [False] * len(lifters)
backtrack(boxes, lifters, 0, steps, used)

# Filter out empty steps
steps = [step for step in steps if step]

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")