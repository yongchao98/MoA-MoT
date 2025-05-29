import heapq

boxes = [256, 53, 286, 205, 261, 209, 242, 143, 354, 55, 72, 161, 268, 305, 43, 246, 146, 331, 84, 128, 88, 358, 71, 93]
lifters = [189, 58, 58, 79, 123, 192, 116]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters, used):
    """ Try to lift a box using available lifters. """
    lifter_indices = []
    remaining_weight = box
    
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            used[i] = True
            if remaining_weight <= 0:
                return lifter_indices
    
    # If we can't lift the box, reset used lifters
    for i in lifter_indices:
        used[i] = False
    return None

def iterative_deepening(boxes, lifters, max_steps):
    for depth in range(1, max_steps + 1):
        steps = []
        used = [False] * len(lifters)
        if backtrack(boxes[:], lifters, 0, steps, used, depth):
            return steps
    return None

def backtrack(boxes, lifters, step, steps, used, max_depth):
    if not boxes:
        return True
    if step >= max_depth:
        return False
    
    current_step = []
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters), max_depth):
            return True
        steps.pop()
        boxes.extend([box for box, _ in current_step])
    
    return False

steps = iterative_deepening(boxes, lifters, 8)

# Validate and print the steps
if steps:
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within the given constraints.")