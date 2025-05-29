from itertools import combinations

# Box weights and lifter capacities
boxes = [49, 186, 267, 243, 352, 74, 160, 115, 138, 301, 250, 145, 294, 232, 144, 293, 287, 358, 267, 266]
lifters = [60, 92, 94, 68, 72, 79]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization dictionary
memo = {}

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            if not any(used[i] for i in indices) and sum(capacities) >= box:
                return list(indices)
    return None

def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        steps.append(step)
        return True
    if len(steps) >= 12:
        return False
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            for idx in lifter_indices:
                used[idx] = True
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step + [(box, lifter_indices)], steps, used):
                memo[state] = True
                return True
            for idx in lifter_indices:
                used[idx] = False
    
    memo[state] = False
    return False

steps = []
used = [False] * len(lifters)
backtrack(boxes, lifters, [], steps, used)

# Print the steps
output = "\n".join(f"Step {i+1}: {step}" for i, step in enumerate(steps))
print(f"<<<{output}>>>")