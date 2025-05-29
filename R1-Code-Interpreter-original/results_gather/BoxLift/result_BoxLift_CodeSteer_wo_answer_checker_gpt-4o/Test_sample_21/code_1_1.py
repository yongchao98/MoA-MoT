from itertools import combinations

boxes = [38, 97, 21, 92, 195, 187, 95, 140, 37, 85, 25, 165]
lifters = [96, 32, 88, 48]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box and not any(used[i] for i in indices):
                return list(indices)
    return None

def backtrack(steps, current_step, boxes, lifters, used):
    if not boxes:
        return steps[:current_step]
    
    if current_step >= 6:
        return None
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            for idx in lifter_indices:
                used[idx] = True
            steps[current_step].append((box, lifter_indices))
            result = backtrack(steps, current_step + 1, boxes[:i] + boxes[i+1:], lifters, used)
            if result is not None:
                return result
            for idx in lifter_indices:
                used[idx] = False
            steps[current_step].pop()
    
    return None

steps = [[] for _ in range(6)]
used = [False] * len(lifters)
solution = backtrack(steps, 0, boxes, lifters, used)

if solution:
    output = "\n".join(f"Step {i+1}: {step}" for i, step in enumerate(solution) if step)
    print(f"<<<{output}>>>")
else:
    print("No solution found within 6 steps.")