boxes = [99, 159, 55, 86, 215, 56, 271, 167, 212, 87, 78, 211, 246, 153, 227, 287]
lifters = [143, 50, 70, 79, 136]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity + lifter <= box:
            total_capacity += lifter
            lifter_indices.append(i)
            if total_capacity == box:
                break
    return lifter_indices if total_capacity >= box else None

def solve_box_lift(boxes, lifters, max_steps):
    steps = []
    while boxes and len(steps) < max_steps:
        step = []
        used = [False] * len(lifters)
        for box in boxes[:]:
            lifter_indices = can_lift(box, lifters, used)
            if lifter_indices is not None:
                step.append((box, lifter_indices))
                boxes.remove(box)
                for index in lifter_indices:
                    used[index] = True
        if step:
            steps.append(step)
        else:
            break

    if not boxes:
        output = []
        for i, step in enumerate(steps):
            output.append(f"Step {i + 1}: {step}")
        print("<<<" + "\n".join(output) + ">>>")
    else:
        print("No solution found within the given steps.")

solve_box_lift(boxes, lifters, 7)