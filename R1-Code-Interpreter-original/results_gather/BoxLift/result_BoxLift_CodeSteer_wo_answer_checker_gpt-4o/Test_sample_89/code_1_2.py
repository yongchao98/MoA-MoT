boxes = [227, 106, 388, 173, 63, 178, 261, 205, 323, 124, 386, 81, 71, 127, 134, 212, 150, 114, 41, 277, 123, 152, 47, 47]
lifters = [166, 106, 122, 181, 102, 119, 147]

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
            if total_capacity >= box:
                return lifter_indices
    return None

def backtrack(boxes, lifters, step, steps, used, memo):
    if not boxes:
        return steps if len(steps) <= 6 else None
    
    if len(steps) >= 6:
        return None
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            for idx in lifter_indices:
                used[idx] = True
            new_steps = steps + [step + [(box, lifter_indices)]]
            result = backtrack(boxes[:i] + boxes[i+1:], lifters, [], new_steps, used, memo)
            if result is not None:
                memo[state] = result
                return result
            for idx in lifter_indices:
                used[idx] = False
    
    memo[state] = None
    return None

used = [False] * len(lifters)
memo = {}
solution = backtrack(boxes, lifters, [], [], used, memo)

# Print the solution
if solution:
    output = []
    for i, step in enumerate(solution):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 6 steps.")