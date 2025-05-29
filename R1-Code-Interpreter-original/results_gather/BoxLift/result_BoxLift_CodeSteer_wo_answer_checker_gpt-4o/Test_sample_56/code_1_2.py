boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
lifters = [78, 98, 51, 75, 43]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity + lifter <= box:
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
            for idx in lifter_indices:
                used[idx] = True
            steps[step].append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step + 1, steps, used):
                return True
            for idx in lifter_indices:
                used[idx] = False
            steps[step].pop()
    
    return False

steps = [[] for _ in range(5)]
used = [False] * len(lifters)
if backtrack(boxes, lifters, 0, steps, used):
    output = ""
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i+1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("No solution found within 5 steps.")