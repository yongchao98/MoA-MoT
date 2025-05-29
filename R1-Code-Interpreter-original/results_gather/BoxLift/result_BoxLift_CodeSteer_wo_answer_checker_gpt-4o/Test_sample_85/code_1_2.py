boxes = [41, 69, 80, 19, 54, 85, 64, 35]
lifters = [45, 45, 62, 74]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters to lift the box
    lifter_indices = []
    remaining_weight = box
    
    for i, capacity in enumerate(lifters):
        if not used[i] and remaining_weight > 0:
            if capacity <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= capacity
                used[i] = True
    
    if remaining_weight <= 0:
        return lifter_indices
    else:
        # Reset used lifters if the combination doesn't work
        for i in lifter_indices:
            used[i] = False
        return None

def backtrack(boxes, lifters, step, max_steps, used, steps, memo):
    if not boxes:
        return step
    
    if step >= max_steps:
        return float('inf')
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    min_steps = float('inf')
    current_step = []
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        min_steps = min(min_steps, backtrack(boxes, lifters, step + 1, max_steps, [False] * len(lifters), steps, memo))
        steps.pop()
    
    memo[state] = min_steps
    return min_steps

max_steps = 4
steps = []
used = [False] * len(lifters)
memo = {}
min_steps = backtrack(boxes, lifters, 0, max_steps, used, steps, memo)

if min_steps <= max_steps:
    output = ""
    for i, step in enumerate(steps):
        output += f"Step {i+1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("It's not possible to lift all boxes within the allowed number of steps.")