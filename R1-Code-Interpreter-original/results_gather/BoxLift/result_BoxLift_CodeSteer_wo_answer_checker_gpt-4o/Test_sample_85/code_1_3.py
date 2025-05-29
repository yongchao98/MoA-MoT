from itertools import combinations

boxes = [41, 69, 80, 19, 54, 85, 64, 35]
lifters = [45, 45, 62, 74]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try all combinations of lifters to lift the box
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo):
                total_capacity = sum(lifters[i] for i in combo)
                if total_capacity >= box:
                    for i in combo:
                        used[i] = True
                    return list(combo)
    return None

def backtrack(boxes, lifters, step, max_steps, used, steps, dp):
    if not boxes:
        return step
    
    if step >= max_steps:
        return float('inf')
    
    state = (tuple(boxes), tuple(used))
    if state in dp:
        return dp[state]
    
    min_steps = float('inf')
    current_step = []
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        min_steps = min(min_steps, backtrack(boxes, lifters, step + 1, max_steps, [False] * len(lifters), steps, dp))
        steps.pop()
    
    dp[state] = min_steps
    return min_steps

max_steps = 4
steps = []
used = [False] * len(lifters)
dp = {}
min_steps = backtrack(boxes, lifters, 0, max_steps, used, steps, dp)

if min_steps <= max_steps:
    output = ""
    for i, step in enumerate(steps):
        output += f"Step {i+1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("It's not possible to lift all boxes within the allowed number of steps.")