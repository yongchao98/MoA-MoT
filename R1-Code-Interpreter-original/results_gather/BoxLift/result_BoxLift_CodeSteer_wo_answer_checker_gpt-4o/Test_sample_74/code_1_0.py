boxes = [166, 144, 53, 213, 51, 156, 197, 311, 177, 358, 172, 134, 179, 145, 91, 188, 352, 294, 292, 88, 97, 394, 123, 294]
lifters = [88, 185, 145, 195, 147, 145]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Precompute possible combinations of lifters that can lift certain weights
def precompute_lifter_combinations(lifters):
    max_weight = sum(lifters)
    dp = [[] for _ in range(max_weight + 1)]
    dp[0] = [[]]  # Base case: weight 0 can be lifted by no lifters

    for i, lifter in enumerate(lifters):
        for weight in range(max_weight, lifter - 1, -1):
            if dp[weight - lifter]:
                for combination in dp[weight - lifter]:
                    dp[weight].append(combination + [i])

    return dp

lifter_combinations = precompute_lifter_combinations(lifters)

def optimistic_estimate(remaining_boxes, lifters):
    # Calculate the minimum number of steps required by considering the heaviest boxes and strongest lifters
    remaining_boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = 0
    i = 0
    while i < len(remaining_boxes):
        weight = 0
        for lifter in lifters:
            if i < len(remaining_boxes) and weight + lifter >= remaining_boxes[i]:
                weight += lifter
                i += 1
        steps += 1
    return steps

def backtrack(boxes, lifters, step, steps, used_boxes, max_steps):
    if step >= max_steps:
        return False
    
    if all(used_boxes):
        return True
    
    remaining_boxes = [box for i, box in enumerate(boxes) if not used_boxes[i]]
    if step + optimistic_estimate(remaining_boxes, lifters) > max_steps:
        return False  # Prune this branch
    
    used_lifters = [False] * len(lifters)
    current_step = []
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Find a valid combination of lifters for this box
        for combination in lifter_combinations[box]:
            if all(not used_lifters[j] for j in combination):
                current_step.append((box, combination))
                for j in combination:
                    used_lifters[j] = True
                used_boxes[i] = True
                break
    
    if current_step:
        steps.append(current_step)
        if backtrack(boxes, lifters, step + 1, steps, used_boxes, max_steps):
            return True
        steps.pop()
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            used_boxes[i] = False
    
    return False

steps = []
used_boxes = [False] * len(boxes)
max_steps = 8

if backtrack(boxes, lifters, 0, steps, used_boxes, max_steps):
    output = []
    for idx, step in enumerate(steps):
        step_str = f"Step {idx + 1}: " + str(step)
        output.append(step_str)
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within the given step limit.")