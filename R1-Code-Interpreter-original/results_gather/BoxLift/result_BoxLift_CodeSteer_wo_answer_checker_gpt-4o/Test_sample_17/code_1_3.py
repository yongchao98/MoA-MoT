from itertools import combinations

boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
lifters = [76, 78, 96, 122, 74, 80]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, selected_lifters = zip(*combo)
            if sum(selected_lifters) >= box and all(not used[i] for i in indices):
                return indices
    return None

def branch_and_bound(boxes, lifters, steps, best_steps, memo):
    if not boxes:
        if len(steps) < len(best_steps[0]):
            best_steps[0] = steps[:]
        return
    
    if len(steps) >= len(best_steps[0]):
        return
    
    state = (tuple(boxes), tuple(steps))
    if state in memo:
        return
    memo.add(state)
    
    used_lifters = [False] * len(lifters)
    step = []
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            for index in lifter_indices:
                used_lifters[index] = True
            boxes.remove(box)
    
    if step:
        steps.append(step)
        branch_and_bound(boxes, lifters, steps, best_steps, memo)
        steps.pop()
        boxes.extend([b[0] for b in step])

best_steps = [[(None, [])] * 8]  # Start with a large number of steps
memo = set()
branch_and_bound(boxes, lifters, [], best_steps, memo)

# Print the best steps
output = "<<<"
for i, step in enumerate(best_steps[0]):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)