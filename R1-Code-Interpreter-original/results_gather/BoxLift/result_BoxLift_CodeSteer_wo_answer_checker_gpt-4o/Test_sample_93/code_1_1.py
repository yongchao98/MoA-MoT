boxes = [64, 70, 386, 351, 113, 77, 314, 333, 266, 399, 193, 44, 181, 200, 238, 175, 370, 118, 337, 134]
lifters = [140, 115, 159, 147, 129, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            total_capacity += lifter
            lifter_indices.append(i)
            used[i] = True
        if total_capacity >= box:
            return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, max_steps):
    if not boxes:
        return True
    if step >= max_steps:
        return False

    used = [False] * len(lifters)
    current_step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if backtrack(remaining_boxes, lifters, step + 1, steps, max_steps):
            return True
        steps.pop()
    
    return False

steps = []
if backtrack(boxes, lifters, 0, steps, 7):
    output = ""
    for i, step in enumerate(steps):
        output += f"Step {i + 1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("No solution found within 7 steps.")