boxes = [399, 277, 278, 282, 209, 309, 281, 53, 325, 239, 363, 363, 308, 244, 392, 323, 255, 79, 280, 336, 395, 318, 203, 364]
lifters = [125, 145, 57, 73, 81, 70]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

# Initialize variables
best_solution = None
best_steps = float('inf')

def can_lift(box, lifters, used):
    remaining_weight = box
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and remaining_weight > 0:
            if lifter <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= lifter
    return remaining_weight <= 0, lifter_indices

def backtrack(step, boxes, lifters, used, current_steps):
    global best_solution, best_steps
    if not boxes:
        if current_steps < best_steps:
            best_steps = current_steps
            best_solution = step[:]
        return

    if current_steps >= best_steps:
        return

    for i, box in enumerate(boxes):
        can_lift_box, lifter_indices = can_lift(box, lifters, used)
        if can_lift_box:
            for index in lifter_indices:
                used[index] = True
            step.append((box, lifter_indices))
            backtrack(step, boxes[:i] + boxes[i+1:], lifters, used, current_steps + 1)
            step.pop()
            for index in lifter_indices:
                used[index] = False

used = [False] * len(lifters)
backtrack([], boxes, lifters, used, 0)

# Print the best solution
output = []
for i, step in enumerate(best_solution):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")