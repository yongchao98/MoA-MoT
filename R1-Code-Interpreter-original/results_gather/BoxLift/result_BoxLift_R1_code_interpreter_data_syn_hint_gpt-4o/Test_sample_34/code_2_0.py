def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                break
    return remaining_weight <= 0, lifter_indices

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        steps.append(step.copy())
        return True
    
    for i, box in enumerate(boxes):
        can_lift_box, lifter_indices = can_lift(box, lifters, used)
        if can_lift_box:
            for idx in lifter_indices:
                used[idx] = True
            step.append((box, lifter_indices))
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step, steps, used):
                return True
            step.pop()
            for idx in lifter_indices:
                used[idx] = False
    return False

boxes = [399, 277, 278, 282, 209, 309, 281, 53, 325, 239, 363, 363, 308, 244, 392, 323, 255, 79, 280, 336, 395, 318, 203, 364]
lifters = [125, 145, 57, 73, 81, 70]
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used = [False] * len(lifters)
backtrack(boxes, lifters, [], steps, used)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")