boxes = [58, 45, 94, 26, 49, 153, 61, 23, 170, 143, 75, 166]
lifters = [57, 61, 104, 98, 70]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i]:
            total_capacity += lifter
            lifter_indices.append(i)
            if total_capacity >= box:
                return lifter_indices
    return None

def backtrack(steps, current_step, boxes, lifters, max_steps):
    if not boxes:
        return steps if len(steps) <= max_steps else None
    if len(steps) >= max_steps:
        return None

    used = [False] * len(lifters)
    step = []
    remaining_boxes = boxes[:]
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                used[index] = True
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

# Print the solution
output = "<<<"
if solution:
    for i, step in enumerate(solution):
        output += f"Step {i+1}: {step}\n"
else:
    output += "No solution found within 4 steps.\n"
output += ">>>"
print(output)