def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            used[i] = True
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
            steps[step].append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step, steps, used):
                return True
            steps[step].pop()
            for idx in lifter_indices:
                used[idx] = False
    
    if step < 6:
        return backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters))
    
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = [[] for _ in range(7)]
    used = [False] * len(lifters)
    
    if backtrack(boxes, lifters, 0, steps, used):
        return steps
    else:
        return "Not possible to lift all boxes in 7 steps or less."

boxes = [362, 202, 192, 149, 357, 335, 209, 193, 334, 375, 374, 167, 276, 70, 282, 201, 242, 158, 236, 212]
lifters = [174, 130, 182, 186, 80, 149, 87]

steps = assign_lifters_to_boxes(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    if step:
        output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)