from itertools import combinations

boxes = [12, 90, 67, 10, 10, 94, 60, 46]
lifters = [56, 77, 72]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo):
                total_capacity = sum(lifters[i] for i in combo)
                if total_capacity >= box:
                    for i in combo:
                        used[i] = True
                    return list(combo)
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        steps.append(step.copy())
        return True

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step, steps, used):
                return True
            step.pop()
            for index in lifter_indices:
                used[index] = False

    return False

steps = []
used = [False] * len(lifters)
backtrack(boxes, lifters, [], steps, used)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")