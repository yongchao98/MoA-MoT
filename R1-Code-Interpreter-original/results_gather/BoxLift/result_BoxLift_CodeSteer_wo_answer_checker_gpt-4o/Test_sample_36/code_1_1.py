from itertools import combinations

boxes = [85, 162, 147, 83, 142, 96, 200, 172, 151, 77, 59, 39]
lifters = [53, 54, 76, 48, 97]

# Sort boxes in descending order
boxes.sort(reverse=True)

best_solution = None
best_steps = float('inf')

def can_lift(box, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box

def backtrack(steps, current_step, remaining_boxes):
    global best_solution, best_steps
    if not remaining_boxes:
        if len(steps) < best_steps:
            best_solution = steps[:]
            best_steps = len(steps)
        return
    if len(steps) >= best_steps or len(steps) >= 6:
        return

    for box in remaining_boxes:
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(current_step[i] == 0 for i in lifter_indices) and can_lift(box, lifter_indices):
                    new_step = current_step[:]
                    for i in lifter_indices:
                        new_step[i] = 1
                    new_steps = steps[:]
                    new_steps.append((box, list(lifter_indices)))
                    new_remaining_boxes = remaining_boxes[:]
                    new_remaining_boxes.remove(box)
                    backtrack(new_steps, new_step, new_remaining_boxes)

backtrack([], [0] * len(lifters), boxes)

# Format the output
output = []
for i, step in enumerate(best_solution):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")