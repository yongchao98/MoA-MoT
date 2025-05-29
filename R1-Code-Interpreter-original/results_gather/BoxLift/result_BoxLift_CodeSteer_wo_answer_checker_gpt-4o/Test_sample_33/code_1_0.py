boxes = [88, 47, 209, 393, 59, 342, 259, 103, 304, 239, 141, 90, 256, 176, 297, 274, 154, 159, 262, 183]
lifters = [173, 188, 50, 119, 177, 86]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, capacity in enumerate(lifters):
        if not used[i] and total_capacity + capacity <= box:
            total_capacity += capacity
            lifter_indices.append(i)
            if total_capacity >= box:
                return lifter_indices
    return None

def branch_and_bound(steps, current_step, boxes, lifters, used):
    if not boxes:
        return True
    if current_step >= 8:
        return False
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            # Mark lifters as used
            for idx in lifter_indices:
                used[idx] = True
            
            # Add to current step
            steps[current_step].append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            
            if branch_and_bound(steps, current_step, remaining_boxes, lifters, used):
                return True
            
            # Backtrack
            for idx in lifter_indices:
                used[idx] = False
            steps[current_step].pop()
    
    # Move to next step
    if steps[current_step]:
        return branch_and_bound(steps, current_step + 1, boxes, lifters, [False] * len(lifters))
    
    return False

steps = [[] for _ in range(8)]
used = [False] * len(lifters)

if branch_and_bound(steps, 0, boxes, lifters, used):
    output = []
    for i, step in enumerate(steps):
        if step:
            output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 8 steps.")