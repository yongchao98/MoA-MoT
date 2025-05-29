boxes = [260, 319, 78, 313, 164, 297, 367, 392, 175, 245, 200, 45, 386, 124, 169, 97, 359, 273, 175, 238]
lifters = [57, 82, 155, 87, 174, 192, 149]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    lifter_indices = []
    total_capacity = 0
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity < box:
            lifter_indices.append(i)
            total_capacity += lifter
        if total_capacity >= box:
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
            # Mark lifters as used
            for index in lifter_indices:
                used[index] = True
            
            # Proceed to the next box
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step, steps, used):
                steps[step].append((box, lifter_indices))
                return True
            
            # Backtrack
            for index in lifter_indices:
                used[index] = False
    
    # Try next step
    return backtrack(boxes, lifters, step + 1, steps, [False] * len(lifters))

steps = [[] for _ in range(7)]
used = [False] * len(lifters)

if backtrack(boxes, lifters, 0, steps, used):
    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within 7 steps.")