def can_lift(box, lifters, used):
    total_capacity = 0
    indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            total_capacity += lifter
            indices.append(i)
    return total_capacity >= box, indices

def backtrack(boxes, lifters, step, steps, used):
    if step >= 6:
        return False
    if not boxes:
        return True
    
    for i, box in enumerate(boxes):
        can_lift_box, indices = can_lift(box, lifters, used)
        if can_lift_box:
            for index in indices:
                used[index] = True
            steps[step].append((box, indices))
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step + 1, steps, used):
                return True
            steps[step].pop()
            for index in indices:
                used[index] = False
    
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = [[] for _ in range(6)]
    used = [False] * len(lifters)
    
    if backtrack(boxes, lifters, 0, steps, used):
        return steps
    else:
        return []

boxes = [368, 220, 131, 287, 153, 175, 153, 217, 262, 115, 365, 394, 291, 222, 142, 126, 189, 82, 159, 367, 356, 168, 188, 355]
lifters = [200, 91, 128, 198, 196, 119, 188]

steps = assign_lifters_to_boxes(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    if step:
        output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)