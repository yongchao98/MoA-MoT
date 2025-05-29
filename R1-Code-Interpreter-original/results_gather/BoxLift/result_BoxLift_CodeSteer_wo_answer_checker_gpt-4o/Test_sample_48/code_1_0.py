boxes = [295, 145, 158, 66, 383, 291, 156, 366, 271, 121, 289, 388, 276, 247, 68, 248, 167, 252, 126, 159]
lifters = [144, 152, 155, 109, 115, 152, 162]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization dictionary
memo = {}

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
    if step >= 6:
        return False
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            # Mark lifters as used
            for idx in lifter_indices:
                used[idx] = True
            steps[step].append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step + 1, steps, used):
                memo[state] = True
                return True
            # Backtrack
            for idx in lifter_indices:
                used[idx] = False
            steps[step].pop()
    
    memo[state] = False
    return False

steps = [[] for _ in range(6)]
used = [False] * len(lifters)

if backtrack(boxes, lifters, 0, steps, used):
    output = ""
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i+1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("No solution found within 6 steps.")