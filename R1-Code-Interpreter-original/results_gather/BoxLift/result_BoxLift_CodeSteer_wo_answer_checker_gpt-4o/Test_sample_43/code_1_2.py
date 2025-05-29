from itertools import combinations

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo) and sum(lifters[i] for i in combo) >= box:
                return combo
    return None

def backtrack(boxes, lifters, step, steps, used, max_steps, memo):
    if step >= max_steps:
        return False
    if not boxes:
        return True
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    for i, box in enumerate(boxes):
        combo = can_lift(box, lifters, used)
        if combo is not None:
            for i in combo:
                used[i] = True
            steps[step].append((box, list(combo)))
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step + 1, steps, used, max_steps, memo):
                memo[state] = True
                return True
            steps[step].pop()
            for i in combo:
                used[i] = False
    
    memo[state] = False
    return False

def assign_lifters_to_boxes(boxes, lifters, max_steps=9):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = [[] for _ in range(max_steps)]
    used = [False] * len(lifters)
    memo = {}
    
    if backtrack(boxes, lifters, 0, steps, used, max_steps, memo):
        return steps
    else:
        return None

boxes = [247, 236, 255, 64, 58, 144, 218, 200, 320, 87, 256, 111, 379, 318, 367, 88, 136, 125, 172, 93]
lifters = [63, 51, 130, 110, 117, 132]

steps = assign_lifters_to_boxes(boxes, lifters)

if steps:
    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within the given step limit.")