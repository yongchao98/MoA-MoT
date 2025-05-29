boxes = [57, 112, 121, 49, 166, 143, 169, 43, 90, 66, 23, 143]
lifters = [103, 49, 95, 97, 98]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    assigned_lifters = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            total_capacity += lifter
            assigned_lifters.append(i)
            used[i] = True
        if total_capacity >= box:
            return assigned_lifters
    return None

def backtrack(steps, current_step, boxes, lifters, max_steps):
    if current_step >= max_steps:
        return False
    if not boxes:
        return True
    
    used = [False] * len(lifters)
    step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        assigned_lifters = can_lift(box, lifters, used)
        if assigned_lifters is not None:
            step.append((box, assigned_lifters))
            remaining_boxes.remove(box)
    
    if step:
        steps.append(step)
        if backtrack(steps, current_step + 1, remaining_boxes, lifters, max_steps):
            return True
        steps.pop()
    
    return False

steps = []
if backtrack(steps, 0, boxes, lifters, 4):
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 4 steps.")