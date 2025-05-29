def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i]:
            total_capacity += lifter
            lifter_indices.append(i)
            if total_capacity >= box:
                return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 5:
        return False
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            # Mark lifters as used
            for index in lifter_indices:
                used[index] = True
            
            # Remove the box and proceed to the next step
            remaining_boxes = boxes[:i] + boxes[i+1:]
            steps[step].append((box, lifter_indices))
            
            if backtrack(remaining_boxes, lifters, step + 1, steps, used):
                return True
            
            # Backtrack
            for index in lifter_indices:
                used[index] = False
            steps[step].pop()
    
    return False

boxes = [18, 93, 23, 70, 72, 22, 96, 96]
lifters = [60, 46, 61]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = [[] for _ in range(5)]
used = [False] * len(lifters)

if backtrack(boxes, lifters, 0, steps, used):
    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i+1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within 5 steps.")