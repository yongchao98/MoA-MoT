boxes = [44, 70, 139, 184, 194, 67, 188, 33, 149, 194, 132, 117]
lifters = [110, 109, 30, 34]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    remaining_weight = box
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
        if remaining_weight <= 0:
            return lifter_indices
    return None

def backtrack(steps, current_step, boxes, lifters, used, memo):
    if not boxes:
        steps.append(current_step[:])
        return True
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            for idx in lifter_indices:
                used[idx] = True
            current_step.append((box, lifter_indices))
            if backtrack(steps, current_step, boxes[:i] + boxes[i+1:], lifters, used, memo):
                memo[state] = True
                return True
            current_step.pop()
            for idx in lifter_indices:
                used[idx] = False
    
    memo[state] = False
    return False

steps = []
used = [False] * len(lifters)
memo = {}
backtrack(steps, [], boxes, lifters, used, memo)

# Format the output
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")