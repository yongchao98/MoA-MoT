def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                return lifter_indices
    
    return None

def backtrack(boxes, lifters, step, steps, used, memo):
    if not boxes:
        return True
    
    if step >= 7:
        return False
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    current_step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            remaining_boxes.remove(box)
            for i in lifter_indices:
                used[i] = True
    
    if current_step:
        steps.append(current_step)
        if backtrack(remaining_boxes, lifters, step + 1, steps, [False] * len(lifters), memo):
            memo[state] = True
            return True
        steps.pop()
    
    memo[state] = False
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    used = [False] * len(lifters)
    memo = {}
    
    if backtrack(boxes, lifters, 0, steps, used, memo):
        return steps
    else:
        return "Cannot lift all boxes in 7 steps or less."

boxes = [71, 72, 167, 64, 339, 45, 243, 53, 52, 280, 164, 322, 362, 75, 158, 142, 150, 353, 182, 282]
lifters = [139, 55, 172, 82, 129, 59, 130]

steps = assign_lifters_to_boxes(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)