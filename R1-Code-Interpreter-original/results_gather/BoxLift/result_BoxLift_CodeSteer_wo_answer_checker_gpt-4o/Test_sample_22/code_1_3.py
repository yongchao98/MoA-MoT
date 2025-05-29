def can_lift(box, lifters, used):
    total_capacity = 0
    indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            total_capacity += lifter
            indices.append(i)
        if total_capacity >= box:
            return indices
    return None

def backtrack(boxes, lifters, used, step, steps, best_steps):
    # Debugging: Print key decision points
    if step >= len(best_steps):
        return
    if not boxes:
        if step < len(best_steps):
            best_steps.clear()
            best_steps.extend(steps)
        return
    for i, box in enumerate(boxes):
        indices = can_lift(box, lifters, used)
        if indices is not None:
            for index in indices:
                used[index] = True
            steps.append((box, indices))
            backtrack(boxes[:i] + boxes[i+1:], lifters, used, step + 1, steps, best_steps)
            steps.pop()
            for index in indices:
                used[index] = False

def box_lift(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    used = [False] * len(lifters)
    best_steps = [None] * 12  # More than the max allowed steps
    backtrack(boxes, lifters, used, 0, [], best_steps)
    return best_steps

# Test with smaller input
boxes = [215, 113, 62, 255, 246, 148, 283, 149, 295, 250, 114, 48, 245, 145, 153, 50]
lifters = [82, 55, 159, 56, 88]

steps = box_lift(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    if step is not None:
        output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)